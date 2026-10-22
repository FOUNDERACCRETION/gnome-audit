"""Register items 30 (zero-gap provenance), 32 (IPW diagnostics),
37 (direct standardisation of the GNoME/MP aggregate inversion),
38 (all matched cells), 40 (MP sub-hull margins)."""
import ast
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261011)
C="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"

s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId","Elements","NSites","Bandgap",C])
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId",C])
s[C]=pd.to_numeric(s[C],errors="coerce"); r[C]=pd.to_numeric(r[C],errors="coerce")
s["bg"]=pd.to_numeric(s.Bandgap,errors="coerce")
s=s.dropna(subset=[C]); r=r.dropna(subset=[C])
R=set(r.MaterialId); s["rec"]=s.MaterialId.isin(R).astype(float)
els=s["Elements"].map(ast.literal_eval).map(set)
s["n_el"]=els.map(len); s["anion"]=els.map(ac)
s["Ln"]=els.map(lambda e:bool(e&LN)).astype(float)
s["d3"]=els.map(lambda e:bool(e&TM3D)).astype(float)
s=s[s.n_el<=4].copy()

print("=== item 30: exact-zero band gaps in the summary file ===")
mm=s[s.rec==1]
g=mm.bg
print("  matched recomputed rows with <=4 el: %d"%len(mm))
print("  exactly zero gaps: %d"%int((g==0).sum()))
raw=mm.loc[g==0,"Bandgap"].astype(str)
print("  distinct raw forms:", dict(raw.value_counts().head(6)))
print("  gaps in (0, 1e-3): %d ;  min strictly positive gap: %.6g"%(((g>0)&(g<1e-3)).sum(), g[g>0].min()))
print("  negative gaps: %d ;  infinite: %d"%(int((g<0).sum()), int(np.isinf(g).sum())))
zc=mm[g==0]
print("  anion mix of zero-gap rows:", dict(zc.anion.value_counts()))
print("  smallest five positive gaps: %s"%np.sort(g[g>0].values)[:5])

# ---------- item 32: IPW diagnostics ----------
F=pd.DataFrame({
  "Ln":s.Ln.values,"d3":s.d3.values,
  "oxide":(s.anion=="oxide").astype(float).values,
  "pnictide":(s.anion=="pnictide").astype(float).values,
  "halide":(s.anion=="halide").astype(float).values,
  "chalc":(s.anion=="chalcogenide").astype(float).values,
  "bcs":(s.anion=="boride/carbide/silicide").astype(float).values,
  "n_el":(s.n_el-s.n_el.mean()).values/s.n_el.std(),
  "lognsites":(np.log(s.NSites)-np.log(s.NSites).mean()).values/np.log(s.NSites).std(),
  "pbe_e":(s[C]-s[C].mean()).values/s[C].std(),
  "pbe_e2":(((s[C]-s[C].mean())/s[C].std())**2).values,
})
X=np.c_[np.ones(len(F)),F.values]; y=s.rec.values
w=np.zeros(X.shape[1])
for it in range(40):
    p=1/(1+np.exp(-(X@w))); W=np.clip(p*(1-p),1e-9,None)
    step=np.linalg.solve((X*W[:,None]).T@X+1e-8*np.eye(len(w)),X.T@(y-p)); w+=step
    if np.max(np.abs(step))<1e-10: break
ps=1/(1+np.exp(-(X@w)))
s["ps"]=ps
print("\n=== item 32: IPW diagnostics ===")
elig=s[s[C]<=0]                      # the flip target: PBE-stable, <=4 elements
print("  propensity fitted on                      %d entries"%len(s))
print("  eligible flip target (PBE-stable, <=4 el) %d"%len(elig))
print("  recomputed among eligible                 %d"%int(elig.rec.sum()))
for lab,sub in (("ALL fitted entries",s),("eligible PBE-stable",elig),
                ("eligible & recomputed",elig[elig.rec==1]),
                ("eligible, not recomputed",elig[elig.rec==0])):
    q=np.percentile(sub.ps,[0,1,5,25,50,75,95,99,100])
    print("  propensity %-26s min %.4f  p1 %.4f  p25 %.3f  med %.3f  p75 %.3f  p99 %.3f  max %.3f"
          %(lab,q[0],q[1],q[3],q[4],q[5],q[7],q[8]))
rec=elig[elig.rec==1].copy()
for trim in (0.0,0.01,0.02,0.05):
    pcl=np.clip(rec.ps.values,trim if trim>0 else None,None)
    wt=1/pcl
    ess=wt.sum()**2/ (wt**2).sum()
    m2=rec.merge(r.rename(columns={C:"scan"}),on="MaterialId")
    fl=(m2.scan>0).astype(float).values
    pcl2=np.clip(m2.ps.values,trim if trim>0 else None,None); wt2=1/pcl2
    print("  trim p>=%.2f : max weight %8.1f   weight sum %10.1f   Kish ESS %8.0f (%.1f%% of n)   IPW rate %.2f%%  naive %.2f%%"
          %(trim,wt2.max(),wt2.sum(),wt2.sum()**2/(wt2**2).sum(),
            100*(wt2.sum()**2/(wt2**2).sum())/len(wt2),
            np.average(fl,weights=wt2)*100, fl.mean()*100))
wt=1/np.clip(rec.ps.values,0.01,None); o=np.argsort(-wt)
print("  share of total weight held by the heaviest 1%% of entries: %.1f%%"
      %(100*wt[o][:max(1,len(wt)//100)].sum()/wt.sum()))
print("  entries with propensity below the 0.01 clip: %d"%int((rec.ps<0.01).sum()))

# ---------- items 37/38: standardisation and all matched cells ----------
print("\n=== items 37/38: GNoME vs Materials Project, all matched cells ===")
rr=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C])
rr[C]=pd.to_numeric(rr[C],errors="coerce"); rr=rr.dropna(subset=[C])
e2=rr["Elements"].map(ast.literal_eval).map(set)
rr["n_elements"]=e2.map(len); rr["anion_class"]=e2.map(ac)
rr["has_Ln"]=e2.map(lambda e:bool(e&LN)); rr["has_3d"]=e2.map(lambda e:bool(e&TM3D))
rr["flip"]=(rr[C]>0).astype(int)
kb=pd.read_csv("results/known_baseline.csv")
gg=rr.groupby(["n_elements","anion_class","has_Ln","has_3d"]).flip.agg(["size","mean"]).reset_index()
gg.columns=["n_elements","anion_class","has_Ln","has_3d","n_G","flip_G"]
j=kb.merge(gg,on=["n_elements","anion_class","has_Ln","has_3d"],how="inner")
j=j[(j.n>=100)&(j.n_G>=100)].copy()
j["delta"]=(j.flip_G-j.flip_known)*100
# cluster-free binomial-ish bootstrap on both sides, independent resampling
def boot(nk,pk,nG,pG,B=100000):
    a=RNG.binomial(nk,pk,B)/nk; b=RNG.binomial(nG,pG,B)/nG
    d=(b-a)*100; return np.percentile(d,[2.5,97.5])
lo=[];hi=[]
for _,x in j.iterrows():
    l,h=boot(int(x.n),x.flip_known,int(x.n_G),x.flip_G); lo.append(l); hi.append(h)
j["lo"]=lo; j["hi"]=hi
j["excl0"]=(j.lo>0)|(j.hi<0)
j=j.sort_values("delta",ascending=False)
print(j[["n_elements","anion_class","has_Ln","has_3d","n_G","flip_G","n","flip_known","delta","lo","hi","excl0"]]
      .to_string(index=False,float_format=lambda v:f"{v:.3f}"))
print("  cells compared: %d   GNoME worse: %d   better: %d   intervals excluding zero: %d"
      %(len(j),int((j.delta>0).sum()),int((j.delta<0).sum()),int(j.excl0.sum())))
j.to_csv("results/known_all_cells.csv",index=False)

print("\n  --- direct standardisation over the matched cells ---")
wG=j.n_G/j.n_G.sum(); wK=j.n/j.n.sum()
print("  crude GNoME rate over matched cells      %.2f%%"%( (j.flip_G*wG).sum()*100))
print("  crude Known rate over matched cells      %.2f%%"%( (j.flip_known*wK).sum()*100))
print("  GNoME standardised to Known cell mix     %.2f%%"%( (j.flip_G*wK).sum()*100))
print("  Known standardised to GNoME cell mix     %.2f%%"%( (j.flip_known*wG).sum()*100))
print("  class mix, GNoME matched cells:")
for a,v in (j.groupby("anion_class").n_G.sum()/j.n_G.sum()).sort_values(ascending=False).items():
    print("     %-24s %.1f%%"%(a,v*100))
print("  class mix, Known matched cells:")
for a,v in (j.groupby("anion_class").n.sum()/j.n.sum()).sort_values(ascending=False).items():
    print("     %-24s %.1f%%"%(a,v*100))
print("\n  --- anion-class mix of the whole r2SCAN release (item 37) ---")
for a,v in (rr.anion_class.value_counts()/len(rr)).items():
    print("     %-24s %.1f%%"%(a,v*100))

# ---------- item 40: MP margins ----------
print("\n=== item 40: Materials Project sub-hull margins ===")
mg=pd.read_csv("results/mp_gga_raw.csv")
print("  MP GGA rows: %d"%len(mg))
eh=pd.to_numeric(mg.energy_above_hull,errors="coerce")
print("  energy_above_hull: min %.6g  median %.6g  max %.6g"%(eh.min(),eh.median(),eh.max()))
print("  exactly zero (floored, i.e. on the hull): %d (%.1f%%)"%(int((eh<=1e-6).sum()),100*(eh<=1e-6).mean()))
print("  -> sub-hull margins are NOT recoverable from this field for stable entries.")

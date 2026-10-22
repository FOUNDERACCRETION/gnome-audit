"""External review of 11 Oct 2026, items 4A and 4B.

4A: the paper reported the eligible propensity range as "<0.001 to 0.800" while also
    reporting an untrimmed maximum weight of 6.6e5, which implies a minimum propensity of
    about 1.5e-6. Report exact minima in scientific notation, per population.
4B: the hydride-dropping sensitivity was fitted with hydride as its own anion dummy, so its
    baseline was not the Table 10 model (which pools hydrides into the reference). Refit
    with the Table 10 design so the only thing that changes is which rows are included.
"""
import ast
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
UEL=set("V Cr Mn Fe Co Ni W Mo".split())
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
def irls(X,y,iters=40):
    b=np.zeros(X.shape[1])
    for _ in range(iters):
        p=1/(1+np.exp(-np.clip(X@b,-30,30))); W=np.clip(p*(1-p),1e-9,None)
        try: step=np.linalg.solve((X*W[:,None]).T@X+1e-10*np.eye(X.shape[1]),X.T@(y-p))
        except np.linalg.LinAlgError: break
        b=b+step
        if np.max(np.abs(step))<1e-11: break
    return b

# ---------------- 4B: hydride sensitivity on the Table 10 design ----------------
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements","NSites",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["has_U"]=els.map(lambda e:bool(e&UEL)).astype(float)
r["has_3d_nonU"]=els.map(lambda e:bool(e&TM3D) and not bool(e&UEL)).astype(float)
r["has_Ln"]=els.map(lambda e:bool(e&LN)).astype(float)
r["flip"]=(r[C]>0).astype(float)
# Table 10 standardisation: computed once on the full released set, held fixed
NEL_M,NEL_S=r.n_el.mean(),r.n_el.std()
LNS_M,LNS_S=np.log(r.NSites).mean(),np.log(r.NSites).std()
DUM=["oxide","pnictide","halide","chalcogenide","boride/carbide/silicide"]
NAMES=["intercept","has_U","has_3d_nonU","has_Ln"]+["is_"+d[:5] for d in DUM]+["n_el_z","lognsites_z"]
def design(d):
    cols=[np.ones(len(d)),d.has_U.values,d.has_3d_nonU.values,d.has_Ln.values]
    for a in DUM: cols.append((d.anion==a).astype(float).values)
    cols.append(((d.n_el-NEL_M)/NEL_S).values)
    cols.append(((np.log(d.NSites)-LNS_M)/LNS_S).values)
    return np.column_stack(cols)
b_all=irls(design(r),r.flip.values)
sub=r[r.anion!="hydride"]
b_noh=irls(design(sub),sub.flip.values)
print("=== 4B. hydride sensitivity, Table 10 design (hydrides pooled into the reference) ===")
print("   n with hydrides %d   n without %d   hydrides %d (%.2f%%)"
      %(len(r),len(sub),len(r)-len(sub),100*(len(r)-len(sub))/len(r)))
print("   %-14s%16s%20s%10s"%("term","with hydrides","hydrides dropped","delta"))
rows=[]
for t in ("has_U","has_3d_nonU","has_Ln"):
    i=NAMES.index(t); v1=np.exp(b_all[i]); v2=np.exp(b_noh[i])
    print("   %-14s%16.3f%20.3f%+10.3f"%(t,v1,v2,v2-v1))
    rows.append(dict(quantity="t10_or_with_hydrides_"+t,value=round(v1,4)))
    rows.append(dict(quantity="t10_or_hydrides_dropped_"+t,value=round(v2,4)))
    rows.append(dict(quantity="t10_hydride_drop_delta_"+t,value=round(v2-v1,4)))
print("   (Table 10 prints has_U 3.448, has_3d_nonU 1.101, has_Ln 1.023 -- the")
print("    'with hydrides' column must match those to be a sensitivity check on it.)")
# every coefficient, so the caption can say "no coefficient moves by more than X"
mx=0.0; mxn=""
for i,t in enumerate(NAMES):
    d=abs(np.exp(b_noh[i])-np.exp(b_all[i]))
    if d>mx: mx,mxn=d,t
print("   largest odds-ratio shift of ANY term: %s, %+.4f"%(mxn,np.exp(b_noh[NAMES.index(mxn)])-np.exp(b_all[NAMES.index(mxn)])))
rows.append(dict(quantity="t10_hydride_drop_max_shift",value=round(mx,4)))

# ---------------- 4A: exact propensity diagnostics ----------------
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId","Elements","NSites",C])
s[C]=pd.to_numeric(s[C],errors="coerce"); s=s.dropna(subset=[C])
R=set(r.MaterialId); s["rec"]=s.MaterialId.isin(R).astype(float)
e2=s["Elements"].map(ast.literal_eval).map(set)
s["n_el"]=e2.map(len); s["anion"]=e2.map(ac)
s["Ln"]=e2.map(lambda e:bool(e&LN)).astype(float)
s["d3"]=e2.map(lambda e:bool(e&TM3D)).astype(float)
s=s[s.n_el<=4].copy()
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
  "pbe_e2":(((s[C]-s[C].mean())/s[C].std())**2).values})
X=np.c_[np.ones(len(F)),F.values]; y=s.rec.values
w=np.zeros(X.shape[1])
for it in range(40):
    p=1/(1+np.exp(-np.clip(X@w,-500,500))); W=np.clip(p*(1-p),1e-9,None)
    step=np.linalg.solve((X*W[:,None]).T@X+1e-8*np.eye(len(w)),X.T@(y-p)); w+=step
    if np.max(np.abs(step))<1e-10: break
s["ps"]=1/(1+np.exp(-np.clip(X@w,-500,500)))
print("\n=== 4A. propensity diagnostics, exact, by population ===")
elig=s[s[C]<=0]
pops=[("all fitted entries",s),("eligible (PBE-stable, <=4 el)",elig),
      ("eligible and recomputed",elig[elig.rec==1]),
      ("eligible, not recomputed",elig[elig.rec==0])]
for lab,d in pops:
    q=np.percentile(d.ps,[1,25,50,75,99])
    print("   %-32s n=%7d  min %.3e  p1 %.4f  p25 %.4f  med %.4f  p75 %.4f  p99 %.4f  max %.4f"
          %(lab,len(d),d.ps.min(),q[0],q[1],q[2],q[3],q[4],d.ps.max()))
rec=elig[elig.rec==1]
print("   implied untrimmed maximum weight 1/min(ps) on 'eligible and recomputed': %.1f"%(1/rec.ps.min()))
print("   entries below the 1e-2 floor, eligible and recomputed: %d"%int((rec.ps<1e-2).sum()))
print("   entries below 1e-4, eligible and recomputed: %d"%int((rec.ps<1e-4).sum()))
rows+= [dict(quantity="ps_min_eligible_recomputed",value=float("%.6g"%rec.ps.min())),
        dict(quantity="ps_min_eligible",value=float("%.6g"%elig.ps.min())),
        dict(quantity="ps_max_eligible",value=round(float(elig.ps.max()),4)),
        dict(quantity="untrimmed_max_weight",value=round(float(1/rec.ps.min()),1)),
        dict(quantity="below_1e-2_floor",value=int((rec.ps<1e-2).sum()))]
pd.DataFrame(rows).to_csv("results/reconcile_checks.csv",index=False)
print("\nwrote results/reconcile_checks.csv")

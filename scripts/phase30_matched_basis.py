"""External review items 2, 4C and 5: the matched comparison against previously catalogued
materials.

(i)  The GNoME side was computed release-based (every r2SCAN row in the cell) while the
     Materials Project side is a confirmed transition (GGA(+U)-stable to r2SCAN-unstable).
     Recompute the GNoME side on the confirmed-transition basis so the two sides use the
     same outcome definition, and report both.
(ii) Report the quantities needed to reconcile the standardisation: totals per side across
     the cells used, oxide entries within those totals, cell weights, UNSTANDARDISED
     aggregate rates within those same cells, and standardised rates under each weighting.
"""
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
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
e2=r["Elements"].map(ast.literal_eval).map(set)
r["n_elements"]=e2.map(len); r["anion_class"]=e2.map(ac)
r["has_Ln"]=e2.map(lambda e:bool(e&LN)); r["has_3d"]=e2.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(int)
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C])
s=s.rename(columns={C:"pbe"}); s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
r=r.merge(s,on="MaterialId",how="left")
conf=r[r.pbe.notna()&(r.pbe<=0)]          # confirmed-transition basis: matched, PBE-stable
print("release-based GNoME rows %d   confirmed-transition rows %d"%(len(r),len(conf)))
print("release-based overall rate %.2f%%   confirmed-transition overall rate %.2f%%"
      %(r.flip.mean()*100, conf.flip.mean()*100))
kb=pd.read_csv("results/known_baseline.csv")
KEY=["n_elements","anion_class","has_Ln","has_3d"]
def cells(d,tag):
    g=d.groupby(KEY).flip.agg(["size","mean"]).reset_index()
    g.columns=KEY+["n_"+tag,"flip_"+tag]; return g
def boot_ci(nk,pk,nG,pG,B=100000):
    a=RNG.binomial(nk,pk,B)/nk; b=RNG.binomial(nG,pG,B)/nG
    d=(b-a)*100; return np.percentile(d,2.5),np.percentile(d,97.5)

for basis,src in (("release-based",r),("confirmed-transition",conf)):
    j=kb.merge(cells(src,"G"),on=KEY,how="inner")
    j=j[(j.n_G>=500)&(j.n>=100)].copy()
    j["delta"]=(j.flip_G-j.flip_known)*100
    lo=[];hi=[]
    for _,x in j.iterrows():
        a,b=boot_ci(int(x.n),x.flip_known,int(x.n_G),x.flip_G); lo.append(a); hi.append(b)
    j["lo"]=lo; j["hi"]=hi; j["excl0"]=(j.lo>0)|(j.hi<0)
    j=j.sort_values("delta",ascending=False)
    wG=j.n_G/j.n_G.sum(); wK=j.n/j.n.sum()
    crudeG=(j.flip_G*wG).sum()*100; crudeK=(j.flip_known*wK).sum()*100
    stdG_K=(j.flip_G*wK).sum()*100; stdK_G=(j.flip_known*wG).sum()*100
    ox=j.anion_class=="oxide"
    print("\n################ %s GNoME side ################"%basis.upper())
    print("  cells meeting n_G>=500 and n_K>=100 : %d  (GNoME worse %d, better %d, intervals excluding zero %d)"
          %(len(j),int((j.delta>0).sum()),int((j.delta<0).sum()),int(j.excl0.sum())))
    print("  TOTAL entries across those cells     GNoME %7d     Known %7d"%(j.n_G.sum(),j.n.sum()))
    print("  oxide entries within those totals    GNoME %7d (%.1f%%)  Known %7d (%.1f%%)"
          %(j.loc[ox,"n_G"].sum(),100*j.loc[ox,"n_G"].sum()/j.n_G.sum(),
            j.loc[ox,"n"].sum(),100*j.loc[ox,"n"].sum()/j.n.sum()))
    print("  UNSTANDARDISED rate within these cells  GNoME %.2f%%   Known %.2f%%   difference %+.2f pp"
          %(crudeG,crudeK,crudeG-crudeK))
    print("  standardised to the KNOWN cell mix      GNoME %.2f%%   Known %.2f%%   difference %+.2f pp"
          %(stdG_K,crudeK,stdG_K-crudeK))
    print("  standardised to the GNoME cell mix      GNoME %.2f%%   Known %.2f%%   difference %+.2f pp"
          %(crudeG,stdK_G,crudeG-stdK_G))
    B=100000; d1=[];d2=[]
    for _ in range(B):
        pg=RNG.binomial(j.n_G.values,j.flip_G.values)/j.n_G.values
        pk=RNG.binomial(j.n.values,j.flip_known.values)/j.n.values
        d1.append(((pg-pk)*wK.values).sum()*100)
        d2.append(((pg-pk)*wG.values).sum()*100)
    d1=np.array(d1); d2=np.array(d2)
    print("  standardised difference at the KNOWN mix  %+.2f pp [%+.2f, %+.2f]"
          %(d1.mean(),np.percentile(d1,2.5),np.percentile(d1,97.5)))
    print("  standardised difference at the GNoME mix  %+.2f pp [%+.2f, %+.2f]"
          %(d2.mean(),np.percentile(d2,2.5),np.percentile(d2,97.5)))
    print("  cell weights, GNoME / Known, by anion class (%):")
    a=(j.groupby('anion_class').n_G.sum()/j.n_G.sum()*100)
    b=(j.groupby('anion_class').n.sum()/j.n.sum()*100)
    for k in a.index: print("     %-24s %5.1f  %5.1f"%(k,a[k],b[k]))
    j.to_csv("results/known_cells_%s.csv"%basis.replace("-","_"),index=False)
print("\n=== anion-class mix, for the aggregate comparison (item 37) ===")
print("  full r2SCAN release (n=%d):"%len(r))
for k,v in (r.anion_class.value_counts()/len(r)*100).items(): print("     %-24s %5.1f%%"%(k,v))
print("  all previously catalogued materials in known_baseline.csv (n=%d):"%kb.n.sum())
for k,v in (kb.groupby('anion_class').n.sum()/kb.n.sum()*100).sort_values(ascending=False).items():
    print("     %-24s %5.1f%%"%(k,v))

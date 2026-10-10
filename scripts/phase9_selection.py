"""Referee #2: model recomputation probability, reweight by inverse propensity.
Referee #5: match GNoME vs MP on hull distance as well as chemistry."""
import ast
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
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

s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId","Elements","NSites",C])
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId",C])
s[C]=pd.to_numeric(s[C],errors="coerce"); r[C]=pd.to_numeric(r[C],errors="coerce")
s=s.dropna(subset=[C]); r=r.dropna(subset=[C])
R=set(r.MaterialId); s["rec"]=s.MaterialId.isin(R).astype(float)
els=s["Elements"].map(ast.literal_eval).map(set)
s["n_el"]=els.map(len); s["anion"]=els.map(ac)
s["Ln"]=els.map(lambda e:bool(e&LN)).astype(float)
s["d3"]=els.map(lambda e:bool(e&TM3D)).astype(float)
s=s[s.n_el<=4].copy()     # quinaries are ~0 coverage; propensity undefined
print(f"<=4-element population {len(s)}, recomputed {int(s.rec.sum())} "
      f"({s.rec.mean()*100:.1f}%)")

# propensity model
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
print(f"propensity model converged ({it+1} iters); "
      f"propensity range {ps.min():.3f}-{ps.max():.3f}")
names=["int"]+list(F.columns)
for nm,c in sorted(zip(names,w),key=lambda t:-abs(t[1]))[:6]:
    print(f"   {nm:12s} coef {c:+.3f}  OR {np.exp(c):.3f}")

# IPW flip rate on recomputed set
m=s[s.rec==1].merge(r.rename(columns={C:"scan"}),on="MaterialId")
idx=s.index[s.rec==1]
m["ipw"]=1/np.clip(ps[s.index.get_indexer(idx)],0.01,None)
m["flip"]=(m.scan>0).astype(float)
naive=m.flip.mean(); ipw=np.average(m.flip,weights=m.ipw)
print(f"\nnaive flip (<=4 el): {naive*100:.2f}%   IPW-reweighted: {ipw*100:.2f}%")

# chemistry gradient under IPW
elsm=m["Elements"].map(ast.literal_eval).map(set)
m["anion"]=elsm.map(ac); m["d3b"]=elsm.map(lambda e:bool(e&TM3D))
g=m[m.n_el==4].groupby(["anion","d3b"]).apply(
    lambda d: pd.Series({"n":len(d),"naive":d.flip.mean(),
                         "ipw":np.average(d.flip,weights=d.ipw)}),include_groups=False)
g=g[g.n>=500].sort_values("naive",ascending=False)
print("\n=== quaternary cells: naive vs IPW ===")
print(g.to_string(float_format=lambda x:f"{x:.4f}"))
g.to_csv("results/ipw_cells.csv")
print(f"\nspread naive {(g.naive.max()-g.naive.min())*100:.1f}pp   "
      f"IPW {(g.ipw.max()-g.ipw.min())*100:.1f}pp")

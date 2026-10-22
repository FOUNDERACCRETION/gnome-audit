"""Defect 2: the regressor labelled "open-shell d (U list)" IS the MP U list, which
contains Mo and W (4d/5d, d0 in their common high oxidation states) and excludes Cu
(d9 as Cu2+). The paper then argues that magnetic open-shell d character, not U-list
membership, is the operative variable -- which that regressor cannot test.

Refit with three competing indicators and compare:
  A  U-list membership        V Cr Mn Fe Co Ni W Mo      (what the paper actually fitted)
  B  3d U-list only           V Cr Mn Fe Co Ni           (drop Mo, W)
  C  open-shell 3d            V Cr Mn Fe Co Ni Cu        (drop Mo/W, add Cu)
If the paper's reading is right, C should fit at least as well as A.
"""
import ast
import numpy as np, pandas as pd
C_="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
SETS={"A U-list (V Cr Mn Fe Co Ni W Mo)":set("V Cr Mn Fe Co Ni W Mo".split()),
      "B 3d U-list (drop Mo,W)":set("V Cr Mn Fe Co Ni".split()),
      "C open-shell 3d (+Cu)":set("V Cr Mn Fe Co Ni Cu".split())}
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
def irls(X,y,it=30):
    b=np.zeros(X.shape[1])
    for _ in range(it):
        p=1/(1+np.exp(-X@b)); W=np.clip(p*(1-p),1e-9,None)
        b=b+np.linalg.solve((X*W[:,None]).T@X+1e-9*np.eye(X.shape[1]),X.T@(y-p))
    return b
def ll(X,b,y):
    p=np.clip(1/(1+np.exp(-X@b)),1e-12,1-1e-12); return (y*np.log(p)+(1-y)*np.log(1-p)).sum()
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C_])
r[C_]=pd.to_numeric(r[C_],errors="coerce"); r=r.dropna(subset=[C_])
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["n_el"]=els.map(len)
y=(r[C_]>0).values.astype(float)
anions=[a for a in sorted(r.anion.unique()) if a!="intermetallic"]
base=[np.ones(len(r)),els.map(lambda e:bool(e&LN)).values.astype(float),
      (r.n_el.values-r.n_el.mean())/r.n_el.std()]
for a in anions: base.append((r.anion==a).values.astype(float))
y_mean=y.mean(); null=(y*np.log(y_mean)+(1-y)*np.log(1-y_mean)).sum()
print(f"{'indicator':36}{'n with':>9}{'OR':>8}{'pseudo-R2':>11}{'logLik':>13}")
out={}
for lab,S in SETS.items():
    tm=els.map(lambda e,s=S: bool(e&s)).values.astype(float)
    oth=els.map(lambda e,s=S: bool(e&set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())) and not bool(e&s)).values.astype(float)
    X=np.column_stack(base+[tm,oth]); b=irls(X,y)
    L=ll(X,b,y); out[lab]=L
    print(f"{lab:36}{int(tm.sum()):9d}{np.exp(b[-2]):8.3f}{1-L/null:11.4f}{L:13.1f}")
print("\nHigher log-likelihood = better fit. The paper's reading predicts C >= A.")
best=max(out,key=out.get); print(f"BEST: {best}")
print(f"  A - C log-likelihood difference: {out['A U-list (V Cr Mn Fe Co Ni W Mo)']-out['C open-shell 3d (+Cu)']:+.1f}")
# and the marginal rates of the disputed elements
print("\nflip rate of the disputed elements:")
for sym in ("Mo","W","Cu","Zn","Sc","Ti"):
    h=els.map(lambda e,s=sym: s in e).values
    print(f"  {sym:3} n={int(h.sum()):6d}  flip {y[h].mean()*100:5.1f}%   (global {y.mean()*100:.1f}%)")

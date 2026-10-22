"""Referee round 4, item 1: can additive per-element energy offsets explain the gradient?

ANALYTIC POINT FIRST. A decomposition energy is evaluated at FIXED composition: the
competing phases that define the hull at a compound's composition contain the same
elements in the same total proportions as the compound itself. An additive shift d_e in
element e's reference energy therefore contributes sum_e x_e d_e to the compound AND the
identical amount to the competing combination, and cancels exactly. Elemental reference
errors move FORMATION energies; they cannot move a decomposition energy. The referee's
mechanism is ruled out by construction, not by data.

EMPIRICAL TEST ANYWAY, because the broader question is fair: is Delta merely an additive
function of composition (bookkeeping-like) or does it depend on how elements are combined
(chemistry)? Fit Delta ~ sum_e c_e x_e with no intercept, then ask how much of the
cell-level gradient survives in the residuals, and whether the FLIP-rate gradient survives
after subtracting the fitted additive part from each material's r2SCAN energy.
"""
import ast, re
import numpy as np, pandas as pd
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

s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C])
s=s.rename(columns={C:"pbe"}); s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).rename(columns={C:"scan"})
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["d3"]=els.map(lambda e:bool(e&TM3D)); r["Ln"]=els.map(lambda e:bool(e&LN))
m=r.merge(s,on="MaterialId",how="inner")
m=m[m.pbe<=0].copy()
m["delta"]=(m.scan-m.pbe)*1000.0            # meV/atom
print(f"PBE-stable, matched: {len(m)}   mean Delta {m.delta.mean():+.2f} meV/atom")

# atomic fractions
def frac(c):
    d={}
    for sym,n in re.findall(r"([A-Z][a-z]?)(\d+)",c): d[sym]=d.get(sym,0)+int(n)
    tot=sum(d.values())
    return {k:v/tot for k,v in d.items()}
F=m.Composition.map(frac)
syms=sorted({k for d in F for k in d})
X=np.zeros((len(m),len(syms)),dtype=np.float32)
pos={s_:i for i,s_ in enumerate(syms)}
for i,d in enumerate(F.values):
    for k,v in d.items(): X[i,pos[k]]=v
y=m.delta.values.astype(np.float64)
coef,_,_,_=np.linalg.lstsq(X.astype(np.float64),y,rcond=None)
pred=X@coef
ss_tot=((y-y.mean())**2).sum(); ss_res=((y-pred)**2).sum()
R2=1-ss_res/ss_tot
print(f"\nAdditive per-element model of Delta ({len(syms)} elements, no intercept):")
print(f"  R^2 = {R2:.3f}   residual s.d. {np.std(y-pred):.1f} meV/atom (raw s.d. {y.std():.1f})")
cf=pd.Series(coef,index=syms).sort_values()
print("  most negative (apparent stabilisation): "+", ".join(f"{k} {v:+.0f}" for k,v in cf.head(5).items()))
print("  most positive (apparent destabilisation): "+", ".join(f"{k} {v:+.0f}" for k,v in cf.tail(5).items()))

m["resid"]=y-pred
print("\nCell-level mean Delta, before and after removing the additive part (meV/atom):")
print(f"{'cell':28}{'n':>7}{'meanD':>9}{'resid':>9}")
rows=[]
for (a,d) in [("oxide",True),("pnictide",True),("oxide",False),("pnictide",False),
              ("intermetallic",True),("intermetallic",False)]:
    sub=m[(m.anion==a)&(m.d3==d)]
    if len(sub)<200: continue
    rows.append(dict(anion=a,d3=d,n=len(sub),mean_delta=sub.delta.mean(),mean_resid=sub.resid.mean()))
    print(f"{a+('+3d' if d else ''):28}{len(sub):7d}{sub.delta.mean():+9.1f}{sub.resid.mean():+9.1f}")
pd.DataFrame(rows).to_csv("results/reference_offsets.csv",index=False)

# decisive: does the FLIP-rate gradient survive removing the additive part?
m["scan_adj"]=m.scan-pred/1000.0
m["flip"]=(m.scan>0).astype(int); m["flip_adj"]=(m.scan_adj>0).astype(int)
print("\nFlip rate before / after subtracting the fitted additive composition term:")
print(f"{'cell':28}{'n':>7}{'flip':>8}{'flip_adj':>10}")
for (a,d) in [("oxide",True),("pnictide",True),("intermetallic",False)]:
    sub=m[(m.anion==a)&(m.d3==d)]
    print(f"{a+('+3d' if d else ''):28}{len(sub):7d}{sub.flip.mean()*100:7.1f}%{sub.flip_adj.mean()*100:9.1f}%")
hi=m[(m.anion=="oxide")&m.d3]; lo=m[(m.anion=="intermetallic")&~m.d3]
print(f"\noxide+3d minus plain intermetallic: raw {(hi.flip.mean()-lo.flip.mean())*100:+.1f} pp, "
      f"after removing additive composition {(hi.flip_adj.mean()-lo.flip_adj.mean())*100:+.1f} pp")
print("\nA gradient that survives in the residuals is not an additive-composition artefact.")

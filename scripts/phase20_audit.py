"""Referee round 5 verification of specific challenged numbers."""
import ast
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()); UEL=set("V Cr Mn Fe Co Ni W Mo".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C,"Bandgap"]).rename(columns={C:"pbe"})
s["pbe"]=pd.to_numeric(s.pbe,errors="coerce"); s["bg"]=pd.to_numeric(s.Bandgap,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).rename(columns={C:"scan"})
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["Ln"]=els.map(lambda e:bool(e&LN)); r["n_el"]=els.map(len)
r["flip"]=(r.scan>0).astype(int)
m=r.merge(s[["MaterialId","bg","pbe"]],on="MaterialId",how="left")

print("=== A. BAND GAP ACCOUNTING (defect 5) ===")
matched=m.pbe.notna()
print(f"  recomputed total            {len(m)}")
print(f"  matched to summary file     {matched.sum()}   unmatched {(~matched).sum()}")
bg=m.bg
print(f"  gap finite                  {np.isfinite(bg).sum()}")
print(f"  gap infinite                {np.isinf(bg).sum()}")
print(f"  gap missing (unmatched/NaN) {bg.isna().sum()}")
print(f"  CHECK finite+inf+nan = {np.isfinite(bg).sum()+np.isinf(bg).sum()+bg.isna().sum()}")

print("\n=== B. LANTHANIDE: marginal vs non-Ln (defect 4) ===")
a=r[r.Ln]; b=r[~r.Ln]
print(f"  Ln      n={len(a):7d} ({len(a)/len(r)*100:.1f}% of set)  flip {a.flip.mean()*100:.2f}%")
print(f"  non-Ln  n={len(b):7d}                  flip {b.flip.mean()*100:.2f}%")
print(f"  difference {(a.flip.mean()-b.flip.mean())*100:+.2f} pp  (vs global {r.flip.mean()*100:.2f}%)")
print("  by anion class, Ln minus non-Ln:")
for an in sorted(r.anion.unique()):
    sub=r[r.anion==an]; x=sub[sub.Ln]; y=sub[~sub.Ln]
    if len(x)<200 or len(y)<200: continue
    print(f"    {an:26}{(x.flip.mean()-y.flip.mean())*100:+7.1f} pp   (n_Ln={len(x)})")

print("\n=== C. HALIDE F vs Cl/Br/I BASELINES (defect 27) ===")
hal=r[r.anion=="halide"].copy(); hasF=els[r.anion=="halide"].map(lambda e:"F" in e).values
for lab,sub in (("F-containing",hal[hasF]),("Cl/Br/I only",hal[~hasF])):
    e_sub=els[sub.index]; u=e_sub.map(lambda e:bool(e&UEL)).values
    print(f"  {lab:14} n={len(sub):6d}  overall {sub.flip.mean()*100:5.1f}%   "
          f"U-list present {sub[u].flip.mean()*100:5.1f}% (n={u.sum()})   "
          f"no U-list {sub[~u].flip.mean()*100:5.1f}%")

print("\n=== D. Fe EXCESS ON THE ODDS SCALE (defect 28) ===")
print(f"  {'anion':26}{'with':>8}{'without':>9}{'pp':>7}{'OR':>8}")
for an in sorted(r.anion.unique()):
    sub=r[r.anion==an]; h=els[r.anion==an].map(lambda e:"Fe" in e).values
    if h.sum()<100: continue
    p1=sub[h].flip.mean(); p0=sub[~h].flip.mean()
    orr=(p1/(1-p1))/(p0/(1-p0))
    print(f"  {an:26}{p1*100:7.1f}%{p0*100:8.1f}%{(p1-p0)*100:+7.1f}{orr:8.2f}")

print("\n=== E. MARGINAL RATES BY ELEMENT COUNT (defect 37) ===")
for k in (2,3,4,5):
    sub=r[r.n_el==k]
    if len(sub): print(f"  n_el={k}: n={len(sub):7d} flip {sub.flip.mean()*100:.2f}%")

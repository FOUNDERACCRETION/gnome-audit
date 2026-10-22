"""Referee round 5: one mutually exclusive population accounting for the band gap, and
censoring bounds defined on the eligible release population rather than on all
recomputed entries."""
import ast
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C,"Bandgap","Elements"])
s=s.rename(columns={C:"pbe"}); s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
s["bg"]=pd.to_numeric(s.Bandgap,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).rename(columns={C:"scan"})
m=r.merge(s[["MaterialId","bg","pbe"]],on="MaterialId",how="left")
print("=== BAND-GAP POPULATION ACCOUNTING (mutually exclusive) ===")
unmatched=m.pbe.isna()
mm=m[~unmatched]
inf_=np.isinf(mm.bg); nan_=mm.bg.isna(); fin=np.isfinite(mm.bg)
zero=(mm.bg==0)&fin
pos=fin&(mm.bg>0)
rows=[("Unmatched recomputed IDs (no summary row)",int(unmatched.sum())),
      ("Matched, source gap missing (NaN)",int(nan_.sum())),
      ("Matched, source gap infinite",int(inf_.sum())),
      ("Matched, gap exactly zero",int(zero.sum())),
      ("Matched, finite positive gap",int(pos.sum()))]
tot=0
for lab,n in rows: print(f"  {lab:46}{n:8d}"); tot+=n
print(f"  {'TOTAL':46}{tot:8d}   (recomputed = {len(m)})")
print(f"  finite gaps available (zero + positive)      {int(fin.sum()):8d}")
print("\n  NOTE: unmatched and NaN are NOT infinite and are not insulator evidence.")

print("\n=== CENSORING BOUNDS ON THE ELIGIBLE RELEASE POPULATION ===")
sel=s.dropna(subset=["pbe"]).copy()
sel["n_el"]=sel["Elements"].map(lambda e: len(ast.literal_eval(e)))
elig=sel[(sel.pbe<=0)&(sel.n_el<=4)]
N=len(elig)
inboth=elig.MaterialId.isin(set(r.MaterialId))
n=int(inboth.sum())
obs=elig[inboth].merge(r,on="MaterialId",how="left")
F=int((obs.scan>0).sum())
print(f"  N eligible release entries (PBE-stable, <=4 elements) {N:8d}")
print(f"  n eligible AND recomputed                             {n:8d}")
print(f"  F observed flips among those                          {F:8d}")
print(f"  observed flip rate among recomputed                   {F/n*100:8.2f}%")
print(f"  p_min = F/N                                           {F/N*100:8.2f}%")
print(f"  p_max = (F + (N-n))/N                                 {(F+(N-n))/N*100:8.2f}%")

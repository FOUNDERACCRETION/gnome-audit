"""Q3: does allocating r2SCAN recomputation by chemistry beat uniform/random?
Honest version: estimate cell flip rates on a FIT half, evaluate catch on a
HELD-OUT half. Ranking by rates fitted on the same data would be circular."""
import ast
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
COL="Decomposition Energy Per Atom MP"
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",COL])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e: bool(e&LN)); r["d3"]=els.map(lambda e: bool(e&TM3D))
r[COL]=pd.to_numeric(r[COL],errors="coerce"); r=r.dropna(subset=[COL])
r["flip"]=(r[COL]>0).astype(int)
key=["n_el","anion","Ln","d3"]

idx=RNG.permutation(len(r)); half=len(r)//2
fit=r.iloc[idx[:half]].copy(); hold=r.iloc[idx[half:]].copy()
rates=fit.groupby(key).flip.mean().rename("rate").reset_index()
hold=hold.merge(rates,on=key,how="left")
hold["rate"]=hold["rate"].fillna(fit.flip.mean())

tot=hold.flip.sum(); n=len(hold)
print(f"held-out n={n}, wrong verdicts={tot} ({tot/n*100:.2f}%)\n")
print(f"{'budget':>7} {'chemistry-ranked':>18} {'random':>10} {'lift':>8}")
rows=[]
order=hold.sort_values("rate",ascending=False)
for b in (0.05,0.10,0.20,0.30,0.50):
    k=int(n*b)
    caught=order.flip.values[:k].sum()/tot
    rand=np.mean([hold.flip.sample(k,random_state=int(s)).sum()/tot for s in range(20)])
    rows.append(dict(budget=b,chem=caught,random=rand,lift_pp=(caught-rand)*100))
    print(f"{b*100:6.0f}% {caught*100:17.1f}% {rand*100:9.1f}% {(caught-rand)*100:7.1f}")
pd.DataFrame(rows).to_csv("results/budget_model.csv",index=False)

print("\n=== who are the 22,498 unmatched IDs? ===")
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId"])
inS=r.MaterialId.isin(set(s.MaterialId))
for lab,sub in [("in summary",r[inS]),("NOT in summary",r[~inS])]:
    print(f"\n{lab}: n={len(sub)}, flip={sub.flip.mean()*100:.2f}%")
    print("  n_elements:", sub.n_el.value_counts().sort_index().to_dict())
    print("  top anion :", sub.anion.value_counts().head(3).to_dict())

"""Phase 3a: selection effect. Which GNoME materials got recomputed at r2SCAN,
and does the naive flip rate need reweighting to the full population?"""
import ast
import numpy as np, pandas as pd
RNG = np.random.default_rng(20261009)
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

def feat(df):
    els=df["Elements"].map(ast.literal_eval).map(set)
    df["n_el"]=els.map(len); df["anion"]=els.map(ac)
    df["Ln"]=els.map(lambda e: bool(e&LN)); df["d3"]=els.map(lambda e: bool(e&TM3D))
    return df

print("loading summary (554k)...", flush=True)
s=pd.read_csv("data/stable_materials_summary.csv", usecols=["MaterialId","Elements",COL])
s=feat(s); s[COL]=pd.to_numeric(s[COL],errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv", usecols=["MaterialId","Elements",COL])
r=feat(r); r[COL]=pd.to_numeric(r[COL],errors="coerce")
r["flip"]=r[COL]>0

s["recomputed"]=s.MaterialId.isin(set(r.MaterialId))
print(f"population={len(s)}  recomputed={s.recomputed.sum()} "
      f"({s.recomputed.mean()*100:.1f}%)", flush=True)

# coverage by cell
cov=(s.groupby(["n_el","anion","Ln","d3"])
       .agg(n_pop=("recomputed","size"), n_rec=("recomputed","sum")).reset_index())
cov["coverage"]=cov.n_rec/cov.n_pop
cov=cov[cov.n_pop>=500].sort_values("coverage")
cov.to_csv("results/selection_coverage.csv",index=False)
print("\n=== coverage by cell (n_pop>=500), lowest 10 and highest 10 ===",flush=True)
print(cov.head(10).to_string(index=False,float_format=lambda x:f"{x:.3f}"),flush=True)
print("...",flush=True)
print(cov.tail(10).to_string(index=False,float_format=lambda x:f"{x:.3f}"),flush=True)

# coverage vs PBE hull distance (is the recomputed set biased toward the hull?)
bins=[-np.inf,-0.2,-0.1,-0.05,-0.02,-0.01,0.0,np.inf]
s["pbe_bin"]=pd.cut(s[COL],bins)
hb=(s.groupby("pbe_bin",observed=True)
      .agg(n_pop=("recomputed","size"),n_rec=("recomputed","sum")).reset_index())
hb["coverage"]=hb.n_rec/hb.n_pop
hb.to_csv("results/selection_by_hull_distance.csv",index=False)
print("\n=== coverage vs PBE decomposition energy ===",flush=True)
print(hb.to_string(index=False,float_format=lambda x:f"{x:.3f}"),flush=True)

# population-reweighted flip rate
key=["n_el","anion","Ln","d3"]
rf=r.merge(s[["MaterialId"]+key].drop_duplicates("MaterialId"),on="MaterialId",
           how="left",suffixes=("","_s"))
cell_flip=rf.groupby(key).agg(flip=("flip","mean"),n=("flip","size")).reset_index()
pop=s.groupby(key).size().rename("w").reset_index()
m=cell_flip.merge(pop,on=key,how="inner")
m=m[m.n>=30]
naive=float(r.flip.mean()); reweighted=float((m.flip*m.w).sum()/m.w.sum())
print(f"\n=== naive flip rate      : {naive*100:.2f}%",flush=True)
print(f"=== reweighted to population: {reweighted*100:.2f}%  "
      f"(cells used: {len(m)}, pop covered: {m.w.sum()})",flush=True)
# worst case: unrecomputed all flip / none flip
nrec=len(r); npop=len(s)
print(f"=== censoring bounds: if all {npop-nrec} unrecomputed flipped -> "
      f"{(r.flip.sum()+(npop-nrec))/npop*100:.1f}% ; if none flipped -> "
      f"{r.flip.sum()/npop*100:.1f}%",flush=True)
m.to_csv("results/reweighted_cells.csv",index=False)

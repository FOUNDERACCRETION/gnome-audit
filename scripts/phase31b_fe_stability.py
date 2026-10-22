"""Review item 3, continued: seed stability of the Fe-in-oxides excess interval and of its
bootstrap standard error, which is the contrast the paper's U argument rests on."""
import ast, re
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r=r[els.map(lambda e:"O" in e)].copy(); els=els[r.index]
r["flip"]=(r[C]>0).astype(float)
r["hasFe"]=els.map(lambda e:"Fe" in e)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
flat=r.flip.values.astype(float); fe=r.hasFe.values
r=r.reset_index(drop=True)
flat=r.flip.values.astype(float); fe=r.hasFe.values
gpos=[np.asarray(v) for v in r.groupby("proto").indices.values()]
print("oxide rows %d, Fe-containing %d, prototypes %d"%(len(r),fe.sum(),len(gpos)),flush=True)
pt=(flat[fe].mean()-flat[~fe].mean())*100
print("point estimate %+.2f pp"%pt,flush=True)
rows=[]
NB=int(__import__("os").environ.get("NB","400"))
for seed in (20261009,1,2,3,4):
    rng=np.random.default_rng(seed); K=len(gpos); out=np.empty(NB)
    for b in range(NB):
        pick=np.concatenate([gpos[i] for i in rng.integers(0,K,K)])
        f=flat[pick]; m=fe[pick]
        out[b]=(f[m].mean()-f[~m].mean())*100
    lo,hi=np.percentile(out,[2.5,97.5]); se=out.std(ddof=1)
    print("seed %-9d [%+.2f, %+.2f]  s.e. %.3f"%(seed,lo,hi,se),flush=True)
    rows.append(dict(seed=seed,lo=round(lo,3),hi=round(hi,3),se=round(se,3)))
d=pd.DataFrame(rows)
print("lower endpoint range %.3f pp, upper endpoint range %.3f pp, s.e. range %.3f"
      %(d.lo.max()-d.lo.min(), d.hi.max()-d.hi.min(), d.se.max()-d.se.min()),flush=True)
d["replicates"]=NB
d.to_csv("results/fe_oxide_seed_stability.csv",index=False)
print("wrote results/fe_oxide_seed_stability.csv")

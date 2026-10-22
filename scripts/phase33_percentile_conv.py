"""External review of 11 Oct 2026, item 2. The paper claimed that seed-to-seed movement in
the lower endpoint of the Fe-oxide excess interval is a property of the sampling
distribution rather than of the simulation. That is wrong as stated: a bootstrap percentile
estimator converges as B grows, however diffuse the tail. This measures the convergence.

Resampling is aggregated to the cluster level first -- per cluster, the Fe/non-Fe counts and
flip sums -- so a replicate is a sum over drawn clusters and B can be taken to 400,000."""
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
r["hasFe"]=els.map(lambda e:"Fe" in e).values
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
r=r.reset_index(drop=True)
codes,uniq=pd.factorize(r.proto); G=len(uniq)
nFe=np.zeros(G); sFe=np.zeros(G); nNo=np.zeros(G); sNo=np.zeros(G)
for g,f,h in zip(codes,r.flip.values,r.hasFe.values):
    if h: nFe[g]+=1; sFe[g]+=f
    else: nNo[g]+=1; sNo[g]+=f
print("oxide rows %d, Fe-containing %d, prototypes %d"%(len(r),int(nFe.sum()),G),flush=True)
pt=(sFe.sum()/nFe.sum()-sNo.sum()/nNo.sum())*100
print("point estimate %+.3f pp"%pt,flush=True)
def endpoints(seed,B):
    rng=np.random.default_rng(seed)
    out=np.empty(B); used=0
    CH=50000
    done=0
    while done<B:
        k=min(CH,B-done)
        pick=rng.integers(0,G,(k,G))
        a=nFe[pick].sum(axis=1); b=sFe[pick].sum(axis=1)
        c=nNo[pick].sum(axis=1); d=sNo[pick].sum(axis=1)
        ok=(a>0)&(c>0)
        v=np.where(ok,(b/np.maximum(a,1)-d/np.maximum(c,1))*100,np.nan)
        v=v[ok]
        out[used:used+len(v)]=v; used+=len(v); done+=k
    out=out[:used]
    return np.percentile(out,2.5), np.percentile(out,97.5), out.std(ddof=1)
SEEDS=(20261009,1,2,3,4)
rows=[]
print("\n%8s %22s %22s %10s"%("B","2.5%ile across seeds","97.5%ile across seeds","spread lo/hi"),flush=True)
for B in (400,4000,40000,100000,400000):
    los=[];his=[]
    for s in SEEDS:
        lo,hi,sd=endpoints(s,B); los.append(lo); his.append(hi)
    print("%8d   [%+.3f, %+.3f]      [%+.3f, %+.3f]     %.3f / %.3f"
          %(B,min(los),max(los),min(his),max(his),max(los)-min(los),max(his)-min(his)),flush=True)
    rows.append(dict(replicates=B,lo_min=round(min(los),3),lo_max=round(max(los),3),
                     hi_min=round(min(his),3),hi_max=round(max(his),3),
                     lo_spread=round(max(los)-min(los),3),hi_spread=round(max(his)-min(his),3),
                     lo_mean=round(float(np.mean(los)),3),hi_mean=round(float(np.mean(his)),3)))
d=pd.DataFrame(rows); d.to_csv("results/percentile_convergence.csv",index=False)
print("\nratio of spreads, successive tenfold increases in B:",flush=True)
for i in range(1,len(rows)):
    print("   B %6d -> %6d : lower spread x%.2f (1/sqrt(10)=0.316 expected)"
          %(rows[i-1]["replicates"],rows[i]["replicates"],
            rows[i]["lo_spread"]/rows[i-1]["lo_spread"]),flush=True)
print("\nwrote results/percentile_convergence.csv")

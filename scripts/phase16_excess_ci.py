"""Cluster-robust intervals on the within-anion excess table (Section 4.3).

The paper reasons from the Fe oxide excess of +5.0 pp against +29 to +47 outside, so that
difference needs an interval, not a point estimate. Resample structural prototypes with
replacement, recomputing the with-minus-without difference inside each anion class.

External review of 11 October 2026 asked whether the printed endpoints are stable to Monte
Carlo variation, and at 400 replicates they were not: the lower endpoint of the Fe oxide
interval moved 0.9 points between seeds. A bootstrap percentile estimator converges as the
replicate count grows, so the fix is replicates, not an explanation. Resampling is therefore
aggregated to the cluster level -- per prototype, the with-element and without-element counts
and flip sums -- which turns a replicate into a sum over drawn clusters and makes 100,000
replicates cheaper than the old implementation's 400.
phase33_percentile_conv.py measures the convergence directly.
"""
import ast, os, re
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))
NB=int(os.environ.get("NB","100000"))
SEED=int(os.environ.get("SEED","20261009"))
CHUNK=int(os.environ.get("CHUNK","4000"))
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).reset_index(drop=True)
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["flip"]=(r[C]>0).astype(float)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
ELS=("Fe","Co","Ni")
for sym in ELS: r["has_"+sym]=els.map(lambda e,s=sym: s in e).values
rows=[]
print(f"{'el':4}{'anion':26}{'n_with':>8}{'excess':>9}{'95% CI':>18}{'clusters':>10}"
      f"   ({NB} replicates, seed {SEED})",flush=True)
for anion in sorted(r.anion.unique()):
    sub=r[r.anion==anion].reset_index(drop=True)
    codes,uniq=pd.factorize(sub.proto); G=len(uniq)
    fl=sub.flip.values
    agg={}
    for sym in ELS:
        h=sub["has_"+sym].values
        nw=np.zeros(G); sw=np.zeros(G); no=np.zeros(G); so=np.zeros(G)
        np.add.at(nw,codes[h],1.0);  np.add.at(sw,codes[h],fl[h])
        np.add.at(no,codes[~h],1.0); np.add.at(so,codes[~h],fl[~h])
        agg[sym]=(nw,sw,no,so,int(h.sum()))
    rng=np.random.default_rng(SEED)
    draws={sym:[] for sym in ELS}
    done=0
    while done<NB:
        k=min(CHUNK,NB-done)
        pick=rng.integers(0,G,(k,G))
        for sym in ELS:
            nw,sw,no,so,cnt=agg[sym]
            if cnt<100: continue
            a=nw[pick].sum(axis=1); b=sw[pick].sum(axis=1)
            c=no[pick].sum(axis=1); d=so[pick].sum(axis=1)
            ok=(a>0)&(c>0)
            draws[sym].append(((b[ok]/a[ok])-(d[ok]/c[ok]))*100)
        done+=k
    for sym in ELS:
        nw,sw,no,so,cnt=agg[sym]
        if cnt<100: continue
        bs=np.concatenate(draws[sym])
        obs=(sw.sum()/nw.sum()-so.sum()/no.sum())*100
        lo,hi=np.percentile(bs,[2.5,97.5])
        rows.append(dict(element=sym,anion=anion,n_with=cnt,excess_pp=round(obs,1),
                         lo=round(lo,1),hi=round(hi,1),n_clusters=G,
                         replicates=len(bs),se_pp=round(bs.std(ddof=1),2)))
        print(f"{sym:4}{anion:26}{cnt:8d}{obs:+9.1f}   [{lo:+6.1f},{hi:+6.1f}]{G:10d}",flush=True)
rows.sort(key=lambda d:(ELS.index(d["element"]),d["anion"]))
pd.DataFrame(rows).to_csv("results/within_anion_excess_ci.csv",index=False)
fe=[x for x in rows if x["element"]=="Fe"]
ox=[x for x in fe if x["anion"]=="oxide"][0]
bc=[x for x in fe if x["anion"]=="boride/carbide/silicide"][0]
print(f"\nKEY CONTRAST: Fe excess in oxides {ox['excess_pp']:+.1f} [{ox['lo']},{ox['hi']}] "
      f"against borides {bc['excess_pp']:+.1f} [{bc['lo']},{bc['hi']}]",flush=True)
print("-> results/within_anion_excess_ci.csv",flush=True)

"""Review item 3, final part: which elements cross the Benjamini-Hochberg threshold, and
whether the same ones do so under different bootstrap seeds. Aggregated to the cluster
level first, so the resampling is a sum over drawn clusters rather than a row gather."""
import ast, re
from math import erf, sqrt
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).reset_index(drop=True)
els=r["Elements"].map(ast.literal_eval).map(set)
r["flip"]=(r[C]>0).astype(float)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
from collections import Counter
cc=Counter(e for s in els for e in s)
KEEP=sorted([e for e,n in cc.items() if n>=300])
print("elements with n>=300: %d"%len(KEEP),flush=True)
ei={e:i for i,e in enumerate(KEEP)}
codes,uniq=pd.factorize(r.proto); G=len(uniq)
CNT=np.zeros((G,len(KEEP))); SUM=np.zeros((G,len(KEEP)))
TOTN=np.zeros(G); TOTS=np.zeros(G)
fl=r.flip.values
for i,(g,s) in enumerate(zip(codes,els)):
    TOTN[g]+=1; TOTS[g]+=fl[i]
    for e in s:
        j=ei.get(e)
        if j is not None: CNT[g,j]+=1; SUM[g,j]+=fl[i]
print("clusters %d"%G,flush=True)
obs=SUM.sum(axis=0)/CNT.sum(axis=0) - TOTS.sum()/TOTN.sum()
import os
NB=int(os.environ.get("NB","20000"))
CH=int(os.environ.get("CH","500"))
def bh(seed,B=None):
    """Resample prototypes with replacement. A draw of G clusters with replacement has the
    same distribution as multinomial counts over the G clusters, and counts turn each
    replicate into one matrix-vector product, so B can be raised far above 300."""
    B=NB if B is None else B
    rng=np.random.default_rng(seed); blocks=[]
    p=np.full(G,1.0/G); done=0
    while done<B:
        k=min(CH,B-done)
        W=rng.multinomial(G,p,size=k).astype(float).T          # (G, k)
        c=CNT.T@W; s=SUM.T@W                                   # (83, k)
        gm=(TOTS@W)/(TOTN@W)                                   # (k,)
        with np.errstate(invalid="ignore",divide="ignore"):
            blocks.append(np.where(c>0,s/np.where(c>0,c,1),np.nan)-gm[None,:])
        done+=k
    D=np.concatenate(blocks,axis=1).T                          # (B, 83)
    se=np.nanstd(D,axis=0,ddof=1)
    z=np.abs(obs)/np.where(se>0,se,np.nan)
    p=np.array([2*(1-0.5*(1+erf(zz/sqrt(2)))) if np.isfinite(zz) else 1.0 for zz in z])
    o=np.argsort(p); m=len(p)
    adj=np.minimum.accumulate((p[o]*m/np.arange(1,m+1))[::-1])[::-1]
    sig=set(np.array(KEEP)[o][adj<0.05])
    padj=dict(zip(np.array(KEEP)[o],adj))
    return sig,padj
SEEDS=(20261009,1,2,3,4)
sets=[];padjs=[]
for s in SEEDS:
    a,b=bh(s); sets.append(a); padjs.append(b)
    print("seed %-9d significant %d  (B=%d)"%(s,len(a),NB),flush=True)
core=set.intersection(*sets); union=set().union(*sets)
print("\nsignificant under ALL five seeds : %d"%len(core))
print("significant under at least one   : %d"%len(union))
print("elements that move between seeds : %s"%sorted(union-core))
for e in sorted(union-core):
    print("   %-3s adjusted p by seed: %s"%(e,[round(pa[e],4) for pa in padjs]))
base=sets[0]
print("\nin the primary seed's set, excess flip rate of the movers:")
for e in sorted(union-core):
    print("   %-3s observed contrast %+.2f pp, in primary set: %s"%(e,obs[ei[e]]*100,e in base))
pd.DataFrame([dict(quantity="bh_core_all_seeds",value=len(core)),
              dict(quantity="bh_union_any_seed",value=len(union)),
              dict(quantity="bh_movers",value=len(union-core))]).to_csv(
              "results/stability_checks.csv",index=False)
print("\nwrote results/stability_checks.csv")

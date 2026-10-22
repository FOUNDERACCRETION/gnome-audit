"""External review of 11 Oct 2026, item 2, applied everywhere rather than in one table.

A bootstrap percentile estimator converges as the replicate count grows. At the counts this
paper originally used -- 600 for cell intervals, 300 for per-element -- the printed endpoints
carried more Monte Carlo noise than the precision they were printed to: measured across five
seeds, per-element endpoints moved 0.3 to 0.9 points and the headline oxide cell's moved 0.5.
This recomputes every cell interval, every per-element interval and the element-versus-global
contrasts at 20,000 replicates, which brings the seed-to-seed movement to about 0.1 points,
below the printed precision.

It is affordable because the resampling is aggregated to the cluster level first. Drawing G
prototypes with replacement has the same distribution as multinomial counts over the G
prototypes, so a replicate becomes one matrix-vector product against per-cluster count and
flip-sum vectors instead of a row gather. Point estimates are unchanged by construction;
only the intervals move.
"""
import ast, os, re
from math import erf, sqrt
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
CA="Decomposition Energy Per Atom All"
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
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))
NB=int(os.environ.get("NB","20000"))
CH=int(os.environ.get("CH","1000"))
SEED=int(os.environ.get("SEED","20261009"))
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).reset_index(drop=True)
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(float)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
print(f"rows {len(r)}, prototypes {r.proto.nunique()}, {NB} replicates, seed {SEED}",flush=True)

def ci_from_counts(n_vec,s_vec,seed,B=NB,ch=CH):
    """Percentile interval for s/n under a cluster bootstrap, via multinomial counts."""
    G=len(n_vec); rng=np.random.default_rng(seed); p=np.full(G,1.0/G); out=[]; done=0
    while done<B:
        k=min(ch,B-done)
        W=rng.multinomial(G,p,size=k).astype(float).T
        c=n_vec@W; t=s_vec@W
        ok=c>0
        out.append(t[ok]/c[ok])
        done+=k
    v=np.concatenate(out)
    return np.percentile(v,2.5), np.percentile(v,97.5), v.std(ddof=1), len(v)

# ---------------- cell intervals ----------------
KEY=["n_el","anion","Ln","d3"]
codes,uniq=pd.factorize(r.proto); G=len(uniq)
rows=[]
grp=r.groupby(KEY)
for key,sub in grp:
    if len(sub)<500: continue
    cc,uu=pd.factorize(sub.proto); g=len(uu)
    n=np.zeros(g); s=np.zeros(g)
    np.add.at(n,cc,1.0); np.add.at(s,cc,sub.flip.values)
    lo,hi,se,B=ci_from_counts(n,s,SEED)
    rows.append(dict(n_el=key[0],anion=key[1],Ln=key[2],d3=key[3],n=len(sub),n_proto=g,
                     rate=sub.flip.mean(),clust_lo=lo,clust_hi=hi,se=se,replicates=B))
cells=pd.DataFrame(rows).sort_values("rate",ascending=False)
# naive binomial interval and the widening ratio, for the appendix table
z975=1.959963984540054
pn=cells.rate.values; nn=cells.n.values
half=z975*np.sqrt(pn*(1-pn)/nn)
cells["naive_lo"]=pn-half; cells["naive_hi"]=pn+half
cells["widening"]=(cells.clust_hi-cells.clust_lo)/(cells.naive_hi-cells.naive_lo)
cells.to_csv("results/cell_ci_converged.csv",index=False)
# the schema figures3.py and the appendix table consume
cells.rename(columns={"Ln":"Ln","d3":"d3"})[
    ["n_el","anion","Ln","d3","n","n_proto","rate","naive_lo","naive_hi",
     "clust_lo","clust_hi","widening"]].to_csv("results/cluster_bootstrap.csv",index=False)
with open("results/appendixA.tex","w") as fh:
    for _,x in cells.iterrows():
        fh.write("%d & %s & %s & %s & %d & %.1f\\%% & [%.1f, %.1f] & %.2f \\\\\n"
                 %(x.n_el,x.anion,"yes" if x.Ln else "no","yes" if x.d3 else "no",
                   x.n,x.rate*100,x.clust_lo*100,x.clust_hi*100,x.widening))
print("wrote results/appendixA.tex and results/cluster_bootstrap.csv",flush=True)
print("cells with n>=500: %d"%len(cells),flush=True)
for _,x in cells.head(3).iterrows():
    print("   %d-el %-24s Ln=%-5s 3d=%-5s  %.4f [%.4f, %.4f]"
          %(x.n_el,x.anion,x.Ln,x.d3,x.rate,x.clust_lo,x.clust_hi),flush=True)
hi_=cells[(cells.n_el==4)&(cells.anion=="oxide")&(~cells.Ln)&(cells.d3)].iloc[0]
lo_=cells[(cells.n_el==4)&(cells.anion=="intermetallic")&(~cells.Ln)&(~cells.d3)].iloc[0]
print("HEADLINE  oxide+3d quaternary %.2f%% [%.2f, %.2f]   plain intermetallic %.2f%% [%.2f, %.2f]"
      %(hi_.rate*100,hi_.clust_lo*100,hi_.clust_hi*100,
        lo_.rate*100,lo_.clust_lo*100,lo_.clust_hi*100),flush=True)

# ---------------- per-element rates, intervals, contrasts, BH ----------------
from collections import Counter
cc=Counter(e for s in els for e in s)
KEEP=sorted([e for e,n in cc.items() if n>=300])
ei={e:i for i,e in enumerate(KEEP)}
CNT=np.zeros((G,len(KEEP))); SUM=np.zeros((G,len(KEEP)))
TOTN=np.zeros(G); TOTS=np.zeros(G)
fl=r.flip.values
for i,(g,s) in enumerate(zip(codes,els)):
    TOTN[g]+=1; TOTS[g]+=fl[i]
    for e in s:
        j=ei.get(e)
        if j is not None: CNT[g,j]+=1; SUM[g,j]+=fl[i]
nproto=(CNT>0).sum(axis=0)
obs=SUM.sum(axis=0)/CNT.sum(axis=0)
gl=TOTS.sum()/TOTN.sum()
rng=np.random.default_rng(SEED); p=np.full(G,1.0/G)
RATE=[]; DIFF=[]; done=0
while done<NB:
    k=min(CH,NB-done)
    W=rng.multinomial(G,p,size=k).astype(float).T
    c=CNT.T@W; s=SUM.T@W; gm=(TOTS@W)/(TOTN@W)
    with np.errstate(invalid="ignore",divide="ignore"):
        rt=np.where(c>0,s/np.where(c>0,c,1),np.nan)
    RATE.append(rt); DIFF.append(rt-gm[None,:])
    done+=k
RATE=np.concatenate(RATE,axis=1); DIFF=np.concatenate(DIFF,axis=1)
lo=np.nanpercentile(RATE,2.5,axis=1); hi=np.nanpercentile(RATE,97.5,axis=1)
se=np.nanstd(DIFF,axis=1,ddof=1)
z=np.abs(obs-gl)/np.where(se>0,se,np.nan)
praw=np.array([2*(1-0.5*(1+erf(zz/sqrt(2)))) if np.isfinite(zz) else 1.0 for zz in z])
o=np.argsort(praw); m=len(praw)
adj=np.empty(m); adj[o]=np.minimum.accumulate((praw[o]*m/np.arange(1,m+1))[::-1])[::-1]
pe=pd.DataFrame(dict(element=KEEP,n=CNT.sum(axis=0).astype(int),n_proto=nproto,
                     flip=obs,clust_lo=lo,clust_hi=hi,se=se,z=z,p_raw=praw,p_BH=adj,
                     sig_BH_05=adj<0.05,replicates=NB)).sort_values("flip",ascending=False)
pe.to_csv("results/per_element_clustered.csv",index=False)
print("\nper-element: %d elements, %d significant after BH at %d replicates"
      %(len(pe),int(pe.sig_BH_05.sum()),NB),flush=True)
for e in ("Fe","Co","Cr","Mn","Hf","Ti","Cu","Ac","Pa","Pm","Tc","Mg","B","Tb"):
    if e in ei:
        x=pe[pe.element==e].iloc[0]
        print("   %-3s n=%6d  %.4f [%.4f, %.4f]  p_BH=%.3g  sig=%s"
              %(e,x.n,x.flip,x.clust_lo,x.clust_hi,x.p_BH,x.sig_BH_05),flush=True)
print("\n-> results/cell_ci_converged.csv, results/per_element_clustered.csv",flush=True)

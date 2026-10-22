"""External review of 11 Oct 2026, item A, and the standing verification request: the
optimised cluster resampling must agree with a direct row-gather resample under equivalent
cluster draws, and a drawn prototype's multiplicity must apply to all of its members,
including across cells.

"Clustering widens every interval" was previously classified as an implementation invariant.
It is not one: cluster-aware uncertainty can be larger or smaller than the independence-based
interval depending on the dependence pattern, the statistic and the resampling design. It is
a property of the recorded run and is now reported as such. What follows is an invariant: two
implementations of the same draw must give the same number.

Exit status 1 if any assertion fails, so run_all treats it as a failed phase.
"""
import ast, re, sys
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
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))
KEY=["n_elements","anion_class","has_Ln","has_3d"]
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).reset_index(drop=True)
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_elements"]=els.map(len); r["anion_class"]=els.map(ac)
r["has_Ln"]=els.map(lambda e:bool(e&LN)); r["has_3d"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(float)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
kb=pd.read_csv("results/known_baseline.csv")
g=r.groupby(KEY).flip.agg(["size","mean"]).reset_index(); g.columns=KEY+["n_G","flip_G"]
j=kb.merge(g,on=KEY,how="inner"); j=j[(j.n_G>=500)&(j.n>=100)].reset_index(drop=True)
NC=len(j)
key_of={tuple(x):i for i,x in enumerate(j[KEY].itertuples(index=False,name=None))}
r["_cell"]=[key_of.get(t,-1) for t in zip(r.n_elements,r.anion_class,r.has_Ln,r.has_3d)]
sub=r[r._cell>=0].reset_index(drop=True)
codes,uniq=pd.factorize(sub.proto); P=len(uniq)
CNT=np.zeros((P,NC)); SUM=np.zeros((P,NC))
np.add.at(CNT,(codes,sub._cell.values),1.0)
np.add.at(SUM,(codes,sub._cell.values),sub.flip.values)
members=[np.where(codes==p)[0] for p in range(P)]
cellv=sub._cell.values; flipv=sub.flip.values
print("cells %d   prototypes %d   materials %d"%(NC,P,len(sub)),flush=True)
fails=[]

print("\n=== 1. aggregate path equals an explicit row-gather resample ===",flush=True)
rng=np.random.default_rng(7)
for trial in range(5):
    W=rng.multinomial(P,np.full(P,1.0/P)).astype(float)
    # aggregate path, as phase37 computes it
    c=CNT.T@W; s=SUM.T@W
    agg=np.where(c>0,s/np.where(c>0,c,1),np.nan)
    # explicit path: materialise the resampled rows, then take per-cell means
    idx=np.concatenate([np.tile(members[p],int(W[p])) for p in range(P) if W[p]>0])
    num=np.zeros(NC); den=np.zeros(NC)
    np.add.at(den,cellv[idx],1.0)
    np.add.at(num,cellv[idx],flipv[idx])
    exp=np.where(den>0,num/np.maximum(den,1),np.nan)
    d=np.nanmax(np.abs(agg-exp))
    ok=(d<1e-12) and np.array_equal(np.isnan(agg),np.isnan(exp))
    print("   trial %d: max |aggregate - explicit| = %.3e   rows %d   %s"
          %(trial,d,len(idx),"ok" if ok else "MISMATCH"),flush=True)
    if not ok: fails.append("aggregate vs explicit, trial %d"%trial)

print("\n=== 2. a drawn prototype's multiplicity applies to all its members, in every cell ===",flush=True)
W=rng.multinomial(P,np.full(P,1.0/P)).astype(float)
idx=np.concatenate([np.tile(members[p],int(W[p])) for p in range(P) if W[p]>0])
obs=np.zeros((P,NC))
np.add.at(obs,(codes[idx],cellv[idx]),1.0)
want=CNT*W[:,None]
d=np.abs(obs-want).max()
print("   max |observed count - multiplicity x base count| = %.3e"%d,flush=True)
if d>1e-9: fails.append("multiplicity not applied consistently")
spanning=np.where((CNT>0).sum(axis=1)>1)[0]
print("   prototypes spanning more than one cell: %d"%len(spanning),flush=True)
bad=0
for p in spanning[:2000]:
    cells=np.where(CNT[p]>0)[0]
    ratios={round(obs[p,cc]/CNT[p,cc],9) for cc in cells}
    if len(ratios)!=1 or abs(list(ratios)[0]-W[p])>1e-9: bad+=1
print("   of those, prototypes whose multiplicity differs between cells: %d"%bad,flush=True)
if bad: fails.append("multiplicity differs between cells for %d prototypes"%bad)

print("\n=== 3. the multinomial draw equals drawing G indices with replacement ===",flush=True)
# distributional check: the mean and variance of a resampled cell rate must agree between
# the two draw mechanisms to within Monte Carlo error at 4,000 replicates each
cell0=int(np.argmax(j.n_G.values))
n0=CNT[:,cell0]; s0=SUM[:,cell0]
B=4000
r1=np.empty(B); r2=np.empty(B)
rg=np.random.default_rng(11)
for b in range(B):
    w=rg.multinomial(P,np.full(P,1.0/P)).astype(float)
    cc=n0@w; r1[b]=(s0@w)/cc if cc>0 else np.nan
for b in range(B):
    pick=rg.integers(0,P,P)
    cc=n0[pick].sum(); r2[b]=s0[pick].sum()/cc if cc>0 else np.nan
m1,m2=np.nanmean(r1),np.nanmean(r2); v1,v2=np.nanstd(r1),np.nanstd(r2)
print("   cell %d: multinomial mean %.5f sd %.5f ; index-draw mean %.5f sd %.5f"
      %(cell0,m1,v1,m2,v2),flush=True)
print("   mean difference %.5f (tol 0.004)   sd ratio %.4f (tol 0.05)"
      %(abs(m1-m2),v1/v2),flush=True)
if abs(m1-m2)>0.004 or abs(v1/v2-1)>0.05:
    fails.append("multinomial and index draws disagree")

print("\n=== 4. recorded-run property, NOT an invariant ===",flush=True)
kc=pd.read_csv("results/known_cells_clustered_release.csv")
w=(kc.hi-kc.lo)/(kc.hi_binom-kc.lo_binom)
print("   in this run the cluster intervals are wider in %d of %d cells (median %.2fx)."
      %(int((w>1).sum()),len(w),w.median()),flush=True)
print("   This is a property of the dependence in these data, not something every",flush=True)
print("   legitimate run must satisfy: a cluster bootstrap can give a narrower interval",flush=True)
print("   than the independence-based one for some statistics and dependence patterns.",flush=True)
pd.DataFrame([dict(quantity="equivalence_checks_failed",value=len(fails)),
              dict(quantity="cells_cluster_wider",value=int((w>1).sum())),
              dict(quantity="cells_compared",value=len(w))]).to_csv(
              "results/resampling_equivalence.csv",index=False)
if fails:
    print("\nFAILED: "+"; ".join(fails),flush=True); sys.exit(1)
print("\nall resampling equivalence checks passed",flush=True)

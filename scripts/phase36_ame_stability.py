"""External review of 11 Oct 2026, item B: the average-marginal-effect intervals were still
at 150 bootstrap refits while everything else went to 20,000-100,000 replicates, and the
response implied a uniform convergence assessment. This measures whether the AME endpoints
are stable to the precision they are printed to, and at what replicate count.

The design matrix is built once and indexed by resampled row positions, rather than rebuilt
from pandas inside the replicate loop, which is what made 150 refits expensive."""
import ast, os, re
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
U8=set("V Cr Mn Fe Co Ni W Mo".split())
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
NB=int(os.environ.get("NB","1000"))
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number","NSites",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).reset_index(drop=True)
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)).astype(float)
r["U"]=els.map(lambda e:bool(e&U8)).astype(float)
r["flip"]=(r[C]>0).astype(float)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
anions=["oxide","pnictide","boride/carbide/silicide","chalcogenide","halide"]
# build the design ONCE
cols=[np.ones(len(r)),r.U.values,r.Ln.values]; nm=["int","has_U","has_Ln"]
for a in anions:
    cols.append((r.anion==a).astype(float).values); nm.append("an_"+a[:6])
cols.append((r.n_el-3.0).values); nm.append("n_el")
cols.append(((np.log(r.NSites)-np.log(r.NSites).mean())/np.log(r.NSites).std()).values); nm.append("lognsites")
AMASK={}
for a in anions:
    m=(r.anion==a).astype(float).values
    cols.append(r.Ln.values*m); nm.append("Lnx"+a[:6]); AMASK["Lnx"+a[:6]]=m
    cols.append(r.U.values*m);  nm.append("Ux"+a[:6]);  AMASK["Ux"+a[:6]]=m
X=np.column_stack(cols); y=r.flip.values
JU=nm.index("has_U"); JL=nm.index("has_Ln")
def irls(Xs,ys,b0=None,iters=25):
    b=np.zeros(Xs.shape[1]) if b0 is None else b0.copy()
    for _ in range(iters):
        p=1/(1+np.exp(-np.clip(Xs@b,-30,30))); W=np.clip(p*(1-p),1e-9,None)
        try: step=np.linalg.solve((Xs*W[:,None]).T@Xs+1e-9*np.eye(Xs.shape[1]),Xs.T@(ys-p))
        except np.linalg.LinAlgError: break
        b=b+step
        if np.max(np.abs(step))<1e-10: break
    return b
def ame(Xs,b,idx,term):
    j=JU if term=="has_U" else JL
    pre="Ux" if term=="has_U" else "Lnx"
    Xa=Xs.copy(); Xb=Xs.copy()
    Xa[:,j]=1.0; Xb[:,j]=0.0
    for k,n_ in enumerate(nm):
        if n_.startswith(pre):
            m=AMASK[n_][idx]
            Xa[:,k]=m; Xb[:,k]=0.0
    pa=1/(1+np.exp(-np.clip(Xa@b,-30,30))); pb=1/(1+np.exp(-np.clip(Xb@b,-30,30)))
    return (pa-pb).mean()*100
b_full=irls(X,y,iters=40)
allidx=np.arange(len(r))
print("point AMEs: has_U %+.3f pp, has_Ln %+.3f pp"
      %(ame(X,b_full,allidx,"has_U"),ame(X,b_full,allidx,"has_Ln")),flush=True)
gpos=[np.asarray(v) for v in r.groupby("proto").indices.values()]
K=len(gpos)
import time
for seed in (int(os.environ.get("SEED","20261011")),):
    rng=np.random.default_rng(seed); t0=time.time()
    du=np.empty(NB); dl=np.empty(NB)
    for b_ in range(NB):
        idx=np.concatenate([gpos[i] for i in rng.integers(0,K,K)])
        Xs=X[idx]; ys=y[idx]
        bb=irls(Xs,ys,b0=b_full,iters=6)
        du[b_]=ame(Xs,bb,idx,"has_U"); dl[b_]=ame(Xs,bb,idx,"has_Ln")
    print("seed %-9d B=%d  has_U [%+.3f, %+.3f]   has_Ln [%+.3f, %+.3f]   %.0fs"
          %(seed,NB,np.percentile(du,2.5),np.percentile(du,97.5),
            np.percentile(dl,2.5),np.percentile(dl,97.5),time.time()-t0),flush=True)

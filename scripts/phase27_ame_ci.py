"""Average marginal effects of the interaction model, with cluster-bootstrap intervals.

External review of 11 Oct 2026 pointed out that these intervals were still at 150 refits
while every other interval in the paper had gone to 20,000-100,000 replicates. At 150 the
endpoints were not stable to the tenth of a point they are printed to. The design matrix is
now built once and indexed by resampled row positions rather than rebuilt from pandas inside
the loop, which makes 1,000 refits affordable; phase36_ame_stability.py measures the
remaining seed-to-seed movement at that count.

Also writes the input checksums (register item 06)."""
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
SEED=int(os.environ.get("SEED","20261011"))
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
pt={"has_U":ame(X,b_full,allidx,"has_U"),"has_Ln":ame(X,b_full,allidx,"has_Ln")}
print("point AMEs: has_U %+.3f pp, has_Ln %+.3f pp"%(pt["has_U"],pt["has_Ln"]),flush=True)
# The paper quotes the per-anion lanthanide AME for halides, oxides and pnictides; point
# estimates are reported for every class but only those three are bootstrapped, which is
# what keeps 1,000 refits inside a sensible runtime.
ANS_PT=anions+["intermetallic"]
ANS=["oxide","pnictide","halide"]
sub_idx={a:np.where((r.anion==a).values)[0] for a in ANS_PT}
pos_of=np.empty(len(r),dtype=np.int64)
def ln_by_anion(Xs,bb,idx,classes=None):
    """Ln AME restricted to each anion class, evaluated on the resampled rows of that class."""
    out={}
    an=r.anion.values[idx]
    for a in (classes or ANS):
        m=(an==a)
        if m.sum()<500: continue
        Xa=Xs[m].copy(); Xb=Xs[m].copy()
        Xa[:,JL]=1.0; Xb[:,JL]=0.0
        for k,n_ in enumerate(nm):
            if n_.startswith("Lnx"):
                mm=AMASK[n_][idx][m]
                Xa[:,k]=mm; Xb[:,k]=0.0
        pa=1/(1+np.exp(-np.clip(Xa@bb,-30,30))); pb=1/(1+np.exp(-np.clip(Xb@bb,-30,30)))
        out[a]=(pa-pb).mean()*100
    return out
pt_anion=ln_by_anion(X,b_full,allidx,ANS_PT)
print("point Ln AME by anion:", {k:round(v,2) for k,v in pt_anion.items()},flush=True)
gpos=[np.asarray(v) for v in r.groupby("proto").indices.values()]
K=len(gpos); rng=np.random.default_rng(SEED)
du=np.empty(NB); dl=np.empty(NB); by={a:[] for a in ANS}
import time; t0=time.time()
for b_ in range(NB):
    idx=np.concatenate([gpos[i] for i in rng.integers(0,K,K)])
    Xs=X[idx]; ys=y[idx]
    bb=irls(Xs,ys,b0=b_full,iters=6)
    du[b_]=ame(Xs,bb,idx,"has_U"); dl[b_]=ame(Xs,bb,idx,"has_Ln")
    for a,v in ln_by_anion(Xs,bb,idx).items(): by[a].append(v)
rows=[dict(term="has_U",ame_pp=round(pt["has_U"],3),lo=round(np.percentile(du,2.5),3),
           hi=round(np.percentile(du,97.5),3),replicates=NB,seed=SEED),
      dict(term="has_Ln",ame_pp=round(pt["has_Ln"],3),lo=round(np.percentile(dl,2.5),3),
           hi=round(np.percentile(dl,97.5),3),replicates=NB,seed=SEED)]
for a,v in by.items():
    if len(v)<10: continue
    v=np.array(v)
    rows.append(dict(term="has_Ln@"+a,ame_pp=round(pt_anion[a],3),
                     lo=round(np.percentile(v,2.5),3),hi=round(np.percentile(v,97.5),3),
                     replicates=NB,seed=SEED))
d=pd.DataFrame(rows); d.to_csv("results/ame_ci.csv",index=False)
print("\n=== cluster-bootstrap intervals on the AMEs (%d refits, seed %d, %.0fs) ==="%(NB,SEED,time.time()-t0),flush=True)
print(d.to_string(index=False),flush=True)
print("\nwrote results/ame_ci.csv",flush=True)
import hashlib
print("\n=== input checksums ===",flush=True)
for f in ("data/stable_materials_summary.csv","data/stable_materials_r2scan.csv"):
    h=hashlib.sha256()
    with open(f,"rb") as fh:
        for chunk in iter(lambda: fh.read(1<<22), b""): h.update(chunk)
    print("  %-40s %s  %d bytes"%(os.path.basename(f),h.hexdigest(),os.path.getsize(f)),flush=True)

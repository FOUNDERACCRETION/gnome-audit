"""Section 3.5's coarser-clustering check, recomputed at the converged replicate count.
The coarse key is anonymised stoichiometry alone, ignoring space group."""
import ast, os, re
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
NB=int(os.environ.get("NB","20000")); CH=1000; SEED=20261009
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).reset_index(drop=True)
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(float)
r["coarse"]=r.Composition.map(anon)
r["primary"]=r["Space Group Number"].astype(str)+"|"+r.coarse
sz=r.groupby("coarse").size()
print("coarse key: %d clusters, mean size %.1f, largest %d"%(len(sz),sz.mean(),sz.max()),flush=True)
print("primary key: %d clusters"%r.primary.nunique(),flush=True)
def ci(sub,key,seed=SEED,B=NB):
    cc,uu=pd.factorize(sub[key]); G=len(uu)
    n=np.zeros(G); s=np.zeros(G)
    np.add.at(n,cc,1.0); np.add.at(s,cc,sub.flip.values)
    rng=np.random.default_rng(seed); p=np.full(G,1.0/G); out=[]; done=0
    while done<B:
        k=min(CH,B-done)
        W=rng.multinomial(G,p,size=k).astype(float).T
        c=n@W; t=s@W; ok=c>0
        out.append(t[ok]/c[ok]); done+=k
    v=np.concatenate(out)
    return np.percentile(v,2.5)*100, np.percentile(v,97.5)*100, G
hi=r[(r.n_el==4)&(r.anion=="oxide")&(~r.Ln)&(r.d3)]
lo=r[(r.n_el==4)&(r.anion=="intermetallic")&(~r.Ln)&(~r.d3)]
print("\n%-34s %8s %22s %22s"%("cell","rate","primary key","coarse key"),flush=True)
for nm,sub in (("4-el oxide, 3d, no Ln",hi),("4-el intermetallic, no Ln/3d",lo)):
    a,b,G1=ci(sub,"primary"); c,d,G2=ci(sub,"coarse")
    print("%-34s %7.2f%%   [%.1f, %.1f] (%d cl)   [%.1f, %.1f] (%d cl)"
          %(nm,sub.flip.mean()*100,a,b,G1,c,d,G2),flush=True)

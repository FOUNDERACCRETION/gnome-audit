"""External review item 3: the seed-stability check covered only the thresholded count of
significant element contrasts. Extend it to the quantities the concern actually names:
bootstrap interval endpoints for the headline cells, the standard error of the key
contrast, and the identity of the elements that cross the significance threshold."""
import ast, re
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
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(float)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
SEEDS=(20261009,1,2,3,4)

def cell_ci(sub,seed,B=600):
    g=[v.values for _,v in sub.groupby("proto").flip]
    K=len(g); rng=np.random.default_rng(seed)
    b=np.array([np.concatenate([g[i] for i in rng.integers(0,K,K)]).mean() for _ in range(B)])
    return np.percentile(b,2.5)*100, np.percentile(b,97.5)*100
hi=r[(r.n_el==4)&(r.anion=="oxide")&(~r.Ln)&(r.d3)]
lo=r[(r.n_el==4)&(r.anion=="intermetallic")&(~r.Ln)&(~r.d3)]
print("=== 1. headline cell interval endpoints across five seeds ===")
print("   cell                                  point      2.5%%ile range    97.5%%ile range")
for nm,sub in (("4-el oxide, 3d, no Ln",hi),("4-el intermetallic, no Ln/3d",lo)):
    los=[];his=[]
    for s in SEEDS:
        a,b=cell_ci(sub,s); los.append(a); his.append(b)
    print("   %-36s %6.2f%%   [%.2f, %.2f]    [%.2f, %.2f]"
          %(nm,sub.flip.mean()*100,min(los),max(los),min(his),max(his)))
    print("        endpoint spread: lower %.3f pp, upper %.3f pp"%(max(los)-min(los),max(his)-min(his)))

print("\n=== 2. Fe excess in oxides: point, interval endpoints and standard error ===")
ox=r[r.anion=="oxide"].copy()
ox["hasFe"]=ox.Elements.map(lambda e:"Fe" in ast.literal_eval(e))
g=[d for _,d in ox.groupby("proto")]
def fe_excess(seed,B=400):
    rng=np.random.default_rng(seed); K=len(g); out=[]
    for _ in range(B):
        d=pd.concat([g[i] for i in rng.integers(0,K,K)])
        a=d[d.hasFe].flip; b=d[~d.hasFe].flip
        if len(a)<1 or len(b)<1: continue
        out.append((a.mean()-b.mean())*100)
    o=np.array(out); return np.percentile(o,2.5),np.percentile(o,97.5),o.std(ddof=1),len(out)
pt=(ox[ox.hasFe].flip.mean()-ox[~ox.hasFe].flip.mean())*100
print("   point estimate %+.2f pp"%pt)
for s in SEEDS:
    a,b,sd,n=fe_excess(s)
    print("   seed %-9d [%+.2f, %+.2f]   bootstrap s.e. %.3f   replicates used %d"%(s,a,b,sd,n))

"""Referee concern #3: are materials independent observations?
GNoME is substitution-generated, so structural near-duplicates may cluster.
Proxy for structural prototype: (space group number, anonymised stoichiometry).
Nd7Os1Pr3Si5 -> counts sorted (1,3,5,7) -> 'A1B3C5D7'."""
import ast, re
from collections import Counter
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
COL="Decomposition Energy Per Atom MP"
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
def anon(comp):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",comp)]
    return "".join(f"{chr(65+i)}{c}" for i,c in enumerate(sorted(cnt)))

r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number",COL])
r[COL]=pd.to_numeric(r[COL],errors="coerce"); r=r.dropna(subset=[COL])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[COL]>0).astype(int)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)

n=len(r); k=r.proto.nunique()
sizes=r.proto.value_counts()
print(f"materials {n}, prototype clusters {k}, mean size {n/k:.2f}")
print(f"largest clusters: {sizes.head(5).to_dict()}")
print(f"singletons {int((sizes==1).sum())} ({(sizes==1).sum()/k*100:.1f}% of clusters)")
print(f"materials in clusters of size>=10: {int(sizes[sizes>=10].sum())} ({sizes[sizes>=10].sum()/n*100:.1f}%)")

# intra-cluster correlation of flip outcome (ANOVA estimator)
g=r.groupby("proto").flip
m=g.mean(); sz=g.size(); pbar=r.flip.mean()
msb=((sz*(m-pbar)**2).sum())/(k-1)
msw=((g.apply(lambda v:((v-v.mean())**2).sum())).sum())/(n-k)
n0=(n-(sz**2).sum()/n)/(k-1)
icc=(msb-msw)/(msb+(n0-1)*msw)
print(f"\nintra-cluster correlation (ICC) of flip outcome: {icc:.4f}")
print(f"mean cluster size n0 = {n0:.2f}  -> design effect 1+(n0-1)*ICC = {1+(n0-1)*icc:.3f}")
print(f"=> CI widths should inflate by about sqrt(deff) = {np.sqrt(1+(n0-1)*icc):.3f}x")

def boot_naive(v,nb=600):
    p=float(np.mean(v)); d=RNG.binomial(len(v),p,nb)/len(v)
    return p,np.percentile(d,2.5),np.percentile(d,97.5)
def boot_clust(sub,nb=600):
    grp=[gg.flip.values for _,gg in sub.groupby("proto")]
    K=len(grp); out=np.empty(nb)
    for b in range(nb):
        idx=RNG.integers(0,K,K)
        out[b]=np.concatenate([grp[i] for i in idx]).mean()
    return float(sub.flip.mean()),np.percentile(out,2.5),np.percentile(out,97.5)

print("\n=== key cells: naive vs cluster bootstrap ===")
print(f"{'cell':<42}{'rate':>7}{'naive CI':>18}{'cluster CI':>20}{'widen':>7}")
rows=[]
for (ne,a,ln,d3),sub in r.groupby(["n_el","anion","Ln","d3"]):
    if len(sub)<500: continue
    p,l1,h1=boot_naive(sub.flip.values); _,l2,h2=boot_clust(sub)
    w=(h2-l2)/(h1-l1)
    rows.append(dict(n_el=ne,anion=a,Ln=ln,d3=d3,n=len(sub),
        n_proto=sub.proto.nunique(),rate=p,naive_lo=l1,naive_hi=h1,
        clust_lo=l2,clust_hi=h2,widening=w))
cb=pd.DataFrame(rows).sort_values("rate",ascending=False)
cb.to_csv("results/cluster_bootstrap.csv",index=False)
for _,x in cb.iterrows():
    lab=f"{x.n_el}el {x.anion[:16]} Ln={str(x.Ln)[0]} 3d={str(x.d3)[0]}"
    print(f"{lab:<42}{x.rate*100:6.1f}%  [{x.naive_lo*100:5.1f},{x.naive_hi*100:5.1f}]"
          f"   [{x.clust_lo*100:5.1f},{x.clust_hi*100:5.1f}]{x.widening:7.2f}")
print(f"\nmedian CI widening factor: {cb.widening.median():.2f}")

# does the headline contrast survive clustering?
top=cb.iloc[0]; bot=cb.iloc[-1]
print(f"\nheadline contrast (n_el matched check below uses 4-element cells only):")
q=cb[cb.n_el==4].sort_values("rate")
lo_cell=q.iloc[0]; hi_cell=q.iloc[-1]
print(f"  highest: {hi_cell.anion} Ln={hi_cell.Ln} 3d={hi_cell.d3} "
      f"{hi_cell.rate*100:.1f}% cluster CI [{hi_cell.clust_lo*100:.1f},{hi_cell.clust_hi*100:.1f}]")
print(f"  lowest : {lo_cell.anion} Ln={lo_cell.Ln} 3d={lo_cell.d3} "
      f"{lo_cell.rate*100:.1f}% cluster CI [{lo_cell.clust_lo*100:.1f},{lo_cell.clust_hi*100:.1f}]")
print(f"  disjoint under clustering: {bool(hi_cell.clust_lo>lo_cell.clust_hi)}")

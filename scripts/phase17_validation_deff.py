"""Two corrections flagged in referee round 4.
(1) Pipeline validation must compare like with like: the source reports 86.8% survival for
    QUATERNARIES and 84% for binaries+ternaries. The paper compared its all-element 13.04%
    against the quaternary figure.
(2) The design effect used the equal-cluster-size formula 1+(n0-1)rho. Cluster sizes here
    run from 1 to 2,953, and with that dispersion the correct version uses the second
    moment: 1 + (sum n_i^2 / sum n_i - 1) * rho.
"""
import ast, re
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["flip"]=(r[C]>0).astype(int)
print("(1) flip rate and implied survival by element count")
for k in (2,3,4,5):
    sub=r[r.n_el==k]
    if len(sub)<50: 
        print(f"   n_el={k}: n={len(sub)} (too few)"); continue
    print(f"   n_el={k}: n={len(sub):7d}  flip {sub.flip.mean()*100:5.2f}%  survival {100-sub.flip.mean()*100:5.2f}%")
bt=r[r.n_el.isin([2,3])]
q=r[r.n_el==4]
print(f"   binaries+ternaries pooled: n={len(bt)}, survival {100-bt.flip.mean()*100:.2f}%  (source reports 84%)")
print(f"   quaternaries:              n={len(q)}, survival {100-q.flip.mean()*100:.2f}%  (source reports 86.8%)")
print(f"   all recomputed:            n={len(r)}, survival {100-r.flip.mean()*100:.2f}%")

def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
g=r.groupby("proto").flip
ni=g.size().values; N=ni.sum(); K=len(ni)
mi=g.mean().values; pbar=r.flip.mean()
# ANOVA ICC
msb=(ni*(mi-pbar)**2).sum()/(K-1)
msw=((r.flip.values-r.proto.map(g.mean()).values)**2).sum()/(N-K)
n0=(N-(ni**2).sum()/N)/(K-1)
rho=(msb-msw)/(msb+(n0-1)*msw)
deff_equal=1+(ni.mean()-1)*rho
deff_unequal=1+((ni**2).sum()/N-1)*rho
print(f"\n(2) clusters K={K}, N={N}, mean size {ni.mean():.2f}, max {ni.max()}")
print(f"   ICC rho = {rho:.4f}")
print(f"   design effect, equal-size formula   1+(nbar-1)rho      = {deff_equal:.2f}  (widening {np.sqrt(deff_equal):.2f}x)")
print(f"   design effect, unequal-size formula 1+(sum n^2/N -1)rho = {deff_unequal:.2f}  (widening {np.sqrt(deff_unequal):.2f}x)")
print(f"   effective n: equal {N/deff_equal:,.0f}   unequal {N/deff_unequal:,.0f}")
print(f"   five largest clusters hold {ni[np.argsort(-ni)][:5].sum()} materials "
      f"({ni[np.argsort(-ni)][:5].sum()/N*100:.1f}% of the data)")

# (3) Does the headline contrast depend on the few enormous clusters?
print("\n(3) headline quaternary contrast, dropping the largest clusters")
order=np.argsort(-ni); keys=g.size().index.values
for drop in (0,5,20,50):
    bad=set(keys[order[:drop]])
    sub=r[~r.proto.isin(bad)]
    e2=sub["Elements"].map(ast.literal_eval).map(set)
    HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
    TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()); LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
    def ac(e):
        if "O" in e: return "oxide"
        if e&HAL: return "halide"
        if e&CHAL: return "chalcogenide"
        if e&PNIC: return "pnictide"
        if "H" in e: return "hydride"
        if e&{"B","C","Si"}: return "boride/carbide/silicide"
        return "intermetallic"
    an=e2.map(ac); d3=e2.map(lambda e:bool(e&TM3D)); ln=e2.map(lambda e:bool(e&LN))
    q4=sub.n_el==4
    hi=sub[q4&(an=="oxide")&d3&~ln]; lo=sub[q4&(an=="intermetallic")&~d3&~ln]
    print(f"   drop {drop:3d} largest ({len(r)-len(sub):6d} materials): "
          f"oxide+3d {hi.flip.mean()*100:5.1f}% (n={len(hi):5d})  "
          f"interm {lo.flip.mean()*100:4.1f}% (n={len(lo):5d})  gap {(hi.flip.mean()-lo.flip.mean())*100:+5.1f}pp")

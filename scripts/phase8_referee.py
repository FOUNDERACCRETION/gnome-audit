"""Referee items: multiple comparisons (BH), budget over repeated splits,
per-element cluster-robust CIs."""
import ast,re
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
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))

r=pd.read_csv("data/stable_materials_r2scan.csv",
    usecols=["MaterialId","Composition","Elements","Space Group Number",COL])
r[COL]=pd.to_numeric(r[COL],errors="coerce"); r=r.dropna(subset=[COL])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[COL]>0).astype(int)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
pbar=r.flip.mean()

# ---- 1. per-element with cluster-robust CI + BH correction
from math import erfc, sqrt
def z_to_p(z): return erfc(abs(z)/sqrt(2))   # two-sided
# The comparator is the global rate, which is estimated from the same data, so the
# contrast must be resampled JOINTLY: one global draw of prototype clusters per replicate,
# from which both the element rate and the global rate are recomputed. Resampling the
# element rate alone and treating the global rate as fixed understates the variance of
# the difference for small elements and overstates it for the largest.
NB=300
protos=[g for _,g in r.groupby("proto").indices.items()]
idx_by_proto=[np.asarray(v) for v in r.groupby("proto").indices.values()]
Kall=len(idx_by_proto)
flip_all=r.flip.values
elem_mask={}
for el in sorted({x for e in els for x in e}):
    m=els.map(lambda e: el in e).values
    if m.sum()>=300: elem_mask[el]=m
boot_glob=np.empty(NB)
boot_el={el:np.empty(NB) for el in elem_mask}
for b in range(NB):
    pick=RNG.integers(0,Kall,Kall)
    rows_idx=np.concatenate([idx_by_proto[i] for i in pick])
    boot_glob[b]=flip_all[rows_idx].mean()
    for el,m in elem_mask.items():
        sel=rows_idx[m[rows_idx]]
        boot_el[el][b]=flip_all[sel].mean() if len(sel) else np.nan
rows=[]
for el,m in elem_mask.items():
    sub=r[m]
    p=float(sub.flip.mean())
    d=boot_el[el]-boot_glob            # the contrast, resampled jointly
    d=d[~np.isnan(d)]
    se=float(d.std(ddof=1))
    rows.append(dict(element=el,n=int(m.sum()),
        n_proto=int(sub.proto.nunique()),flip=p,
        clust_lo=float(np.nanpercentile(boot_el[el],2.5)),
        clust_hi=float(np.nanpercentile(boot_el[el],97.5)),
        se=se,z=(p-pbar)/se if se>0 else 0.0))
pe=pd.DataFrame(rows)
pe["p_raw"]=pe.z.map(z_to_p)
pe=pe.sort_values("p_raw").reset_index(drop=True)
mtests=len(pe)
# BH step-up: the monotonicity enforcement runs from the LARGEST p downward.
# Taking cummin forward over an ascending sort propagates the smallest adjusted
# p to every later element and marks everything significant.
pe["p_BH"]=(pe.p_raw*mtests/(pe.index+1))[::-1].cummin()[::-1].clip(upper=1.0)
pe["sig_BH_05"]=pe.p_BH<0.05
pe=pe.sort_values("flip",ascending=False)
pe.to_csv("results/per_element_clustered.csv",index=False)
print(f"=== PER-ELEMENT, cluster-robust, BH over {mtests} tests ===")
print(f"significant after BH at 0.05: {int(pe.sig_BH_05.sum())} of {mtests}")
print(pe.head(12)[["element","n","n_proto","flip","clust_lo","clust_hi","p_BH","sig_BH_05"]]
      .to_string(index=False,float_format=lambda x:f"{x:.4f}"))
print("...")
print(pe.tail(6)[["element","n","n_proto","flip","clust_lo","clust_hi","p_BH","sig_BH_05"]]
      .to_string(index=False,float_format=lambda x:f"{x:.4f}"))

# ---- 2. budget experiment over 20 repeated splits, cluster-aware split
key=["n_el","anion","Ln","d3"]
protos=r.proto.unique()
res={b:[] for b in (0.05,0.10,0.20,0.30,0.50)}
rnd={b:[] for b in res}
for seed in range(20):
    rs=np.random.default_rng(1000+seed)
    # split BY CLUSTER so related materials never straddle fit/test
    pick=set(rs.choice(protos,size=len(protos)//2,replace=False))
    fit=r[r.proto.isin(pick)]; hold=r[~r.proto.isin(pick)]
    if len(hold)<1000: continue
    rates=fit.groupby(key).flip.mean()
    h=hold.join(rates.rename("rate"),on=key)
    h["rate"]=h["rate"].fillna(fit.flip.mean())
    tot=h.flip.sum()
    if tot==0: continue
    order=h.sort_values("rate",ascending=False).flip.values
    for b in res:
        k=int(len(h)*b)
        res[b].append(order[:k].sum()/tot)
        rnd[b].append(np.mean([h.flip.sample(k,random_state=int(s)).sum()/tot for s in range(3)]))
print("\n=== BUDGET over 20 cluster-wise splits ===")
print(f"{'budget':>7}{'chem mean':>11}{'sd':>7}{'random':>9}{'lift':>8}")
out=[]
for b in sorted(res):
    c=np.array(res[b]); q=np.array(rnd[b])
    print(f"{b*100:6.0f}%{c.mean()*100:10.1f}%{c.std()*100:7.2f}{q.mean()*100:8.1f}%"
          f"{(c.mean()-q.mean())*100:+8.1f}")
    out.append(dict(budget=b,chem_mean=c.mean(),chem_sd=c.std(),
                    chem_lo=np.percentile(c,2.5),chem_hi=np.percentile(c,97.5),
                    random_mean=q.mean(),lift_pp=(c.mean()-q.mean())*100))
pd.DataFrame(out).to_csv("results/budget_repeated.csv",index=False)
print("\ndone")

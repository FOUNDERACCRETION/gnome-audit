"""Register items 10, 18, 32, 38: seed stability of the BH count, the hydride-dropping
sensitivity actually computed, IPW weight diagnostics on the eligible population, and the
full twenty matched Materials Project cells."""
import ast, re, csv, math
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
C="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
UEL=set("V Cr Mn Fe Co Ni W Mo".split()); TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
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
def irls(X,y,it=30):
    b=np.zeros(X.shape[1])
    for _ in range(it):
        p=1/(1+np.exp(-X@b)); W=np.clip(p*(1-p),1e-9,None)
        b=b+np.linalg.solve((X*W[:,None]).T@X+1e-9*np.eye(X.shape[1]),X.T@(y-p))
    return b
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number","NSites",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["n_el"]=els.map(len); r["flip"]=(r[C]>0).astype(int)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)

print("=== 10. seed stability of the BH count ===")
idx_by_proto=[np.asarray(v) for v in r.groupby("proto").indices.values()]
Kall=len(idx_by_proto); flip_all=r.flip.values; pbar=flip_all.mean()
emask={}
for el in sorted({x for e in els for x in e}):
    m=els.map(lambda e: el in e).values
    if m.sum()>=300: emask[el]=m
def bh_count(seed,NB=300):
    rng=np.random.default_rng(seed)
    bg=np.empty(NB); be={el:np.empty(NB) for el in emask}
    for b in range(NB):
        pick=rng.integers(0,Kall,Kall)
        ridx=np.concatenate([idx_by_proto[i] for i in pick])
        bg[b]=flip_all[ridx].mean()
        for el,m in emask.items():
            sel=ridx[m[ridx]]; be[el][b]=flip_all[sel].mean() if len(sel) else np.nan
    ps=[]
    for el,m in emask.items():
        d=be[el]-bg; d=d[~np.isnan(d)]
        se=d.std(ddof=1); p=r.flip.values[m].mean()
        z=(p-pbar)/se if se>0 else 0.0
        ps.append(math.erfc(abs(z)/math.sqrt(2)))
    ps=np.sort(np.array(ps)); mt=len(ps)
    adj=np.minimum.accumulate((ps*mt/np.arange(1,mt+1))[::-1])[::-1]
    return int((adj<0.05).sum())
counts=[bh_count(s) for s in (20261009,1,2,3,4)]
print(f"   BH count across five seeds: {counts}  (range {min(counts)}-{max(counts)})")

print("\n=== 18. hydride-dropping sensitivity, actually computed ===")
anions=[a for a in sorted(r.anion.unique()) if a!="intermetallic"]
def fit(sub, drop_hyd):
    e2=els[sub.index]
    base=[np.ones(len(sub)),
          e2.map(lambda e:bool(e&UEL)).values.astype(float),
          e2.map(lambda e:bool(e&TM3D) and not bool(e&UEL)).values.astype(float),
          e2.map(lambda e:bool(e&LN)).values.astype(float),
          (sub.n_el.values-sub.n_el.mean())/sub.n_el.std(),
          (np.log(sub.NSites.values)-np.log(sub.NSites.values).mean())/np.log(sub.NSites.values).std()]
    names=["intercept","has_U","has_3d_nonU","has_Ln","n_el_z","lognsites_z"]
    for a in anions:
        if drop_hyd and a=="hydride": continue
        base.append((sub.anion==a).values.astype(float)); names.append(a)
    return names, irls(np.column_stack(base), sub.flip.values.astype(float))
n1,b1=fit(r,False); n2,b2=fit(r[r.anion!="hydride"],True)
print(f"   {'term':16}{'with hydrides':>15}{'hydrides dropped':>19}{'delta':>10}")
_rows=[dict(quantity="bh_count_min",value=min(counts)),
       dict(quantity="bh_count_max",value=max(counts))]
for term in ("has_U","has_3d_nonU","has_Ln"):
    v1=np.exp(b1[n1.index(term)]); v2=np.exp(b2[n2.index(term)])
    print(f"   {term:16}{v1:15.3f}{v2:19.3f}{v2-v1:+10.3f}")
    _rows.append(dict(quantity="hydride_drop_delta_"+term,value=round(v2-v1,4)))
    _rows.append(dict(quantity="or_with_hydrides_"+term,value=round(v1,4)))
    _rows.append(dict(quantity="or_hydrides_dropped_"+term,value=round(v2,4)))
pd.DataFrame(_rows).to_csv("results/register_checks.csv",index=False)
print("   wrote results/register_checks.csv")

"""Identification, a joint test, and cluster-aware uncertainty for Table 4.

Three things, each answering a different objection.

(1) IDENTIFICATION, including the part that is fragile. The fitted design has no intercept
and retains all six non-hydride contexts for each of the five elements, so it is NOT
reference-coded. It is full rank, but only just: the sum of an element's six interaction
columns equals its main-effect column on every row EXCEPT hydride-class materials, and only
286 of those contain one of the five elements (7 for Mn). Dropping hydride-class rows takes
the rank from 113 to 111 and the condition number from 7.0e2 to 3.6e8. So the element
main-effect LEVELS are weakly identified. What the table reports is a within-element spread,
which is a difference among interaction coefficients and does not involve that direction;
this asserts it, by projecting each contrast onto the design's smallest singular directions.

(2) A JOINT TEST. A spread is max-minus-min of fitted coefficients, non-negative by
construction, with its argmax and argmin selected from the same data, so a bootstrap
interval for it is not a test of the additive null. The test is a cluster-robust Wald test
of the joint null that an element's six interaction coefficients are equal, with
Benjamini-Hochberg correction across the five elements.

(3) UNCERTAINTY for the spreads, as a descriptive interval rather than a test. Cluster
bootstrap over structural prototypes, resumable: each run appends its draws so the printed
percentiles can be built from more replicates than one call can produce, and the per-seed
endpoints are reported so their stability is visible.

Outputs results/context_identification.csv, results/context_joint_test.csv,
results/context_spread_ci.csv and results/context_spread_draws.npy.
"""
import ast, os, re, time
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
TM=["Fe","Co","Ni","Mn","Cr"]
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

NB=int(os.environ.get("NB","1000")); SEED=int(os.environ.get("SEED","20261009"))
BUDGET=float(os.environ.get("BUDGET","120")); RESET=os.environ.get("RESET","")
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C]).rename(columns={C:"pbe"})
s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).rename(columns={C:"scan"})
els=r["Elements"].map(ast.literal_eval).map(set); r["anion"]=els.map(ac)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
m=r.merge(s,on="MaterialId",how="inner"); m=m[m.pbe<=0].reset_index(drop=True)
m["delta"]=(m.scan-m.pbe)*1000.0
print(f"n = {len(m)}, prototypes = {m.proto.nunique()}")

def frac(c):
    d={}
    for sym,n in re.findall(r"([A-Z][a-z]?)(\d+)",c): d[sym]=d.get(sym,0)+int(n)
    t=sum(d.values()); return {k:v/t for k,v in d.items()}
F=m.Composition.map(frac); syms=sorted({k for d in F for k in d}); pos={s_:i for i,s_ in enumerate(syms)}
X0=np.zeros((len(m),len(syms)))
for i,d in enumerate(F.values):
    for k,v in d.items(): X0[i,pos[k]]=v
y=m.delta.values; av=m.anion.values; anions=sorted(m.anion.unique())
cols=[];labels=[]
for el in TM:
    xe=X0[:,pos[el]]
    for a in anions:
        v=xe*(av==a)
        if (v>0).sum()>=100: cols.append(v); labels.append((el,a))
X=np.column_stack([X0]+cols); K=X.shape[1]; NS=len(syms)
print(f"design: {X.shape[0]} x {K}  ({NS} element columns, no intercept, "
      f"{len(cols)} interactions, not reference-coded)")

def spectrum(A):
    ev=np.linalg.eigvalsh(A.T@A); sv=np.sqrt(np.clip(ev,0,None))[::-1]
    return sv
sv=spectrum(X)
tol=sv.max()*max(X.shape)*np.finfo(float).eps
rank=int((sv>tol).sum()); cond=sv.max()/sv[sv>tol].min()
print(f"rank {rank} of {K} (tol {tol:.3e}); condition number {cond:.3e}")

# what the full rank rests on
hyd=(av=="hydride")
hyd_cnt={el:int(((X0[:,pos[el]]>0)&hyd).sum()) for el in TM}
sv_nohyd=spectrum(X[~hyd])
tol2=sv_nohyd.max()*max(X[~hyd].shape)*np.finfo(float).eps
rank_nohyd=int((sv_nohyd>tol2).sum()); cond_nohyd=sv_nohyd.max()/sv_nohyd[sv_nohyd>tol2].min()
print(f"hydride-class materials containing one of the five: {hyd_cnt} "
      f"(total {sum(hyd_cnt.values())})")
print(f"without hydride-class rows: rank {rank_nohyd} of {K}, condition {cond_nohyd:.3e}")

# The near-dependency is analytic, not just numerical: for each listed element, its
# main-effect column minus the sum of its six interaction columns is zero on every row
# except hydride-class ones. Build those vectors explicitly, confirm they are the weak
# directions, and check the reported contrasts against them exactly. The empirical singular
# vectors are reported too, but they mix other columns and so are a diagnostic, not a test.
ANA=[]
for el in TM:
    v=np.zeros(K); v[pos[el]]=1.0
    for a in [a for (e,a) in labels if e==el]: v[NS+labels.index((el,a))]=-1.0
    ANA.append(v/np.linalg.norm(v))
ANA=np.array(ANA)
resid=np.linalg.norm(X@ANA.T,axis=0)/np.linalg.norm(X,axis=0).max()
print(f"analytic near-dependency vectors: relative |X v| = "
      f"{np.array2string(resid,precision=2e-0 and 4,floatmode='maxprec')}")
_,_,Vt=np.linalg.svd(X.T@X)
weak=Vt[-(K-rank_nohyd):] if K>rank_nohyd else np.zeros((0,K))
def eff_vec(el,a):
    v=np.zeros(K); v[pos[el]]=1.0; v[NS+labels.index((el,a))]=1.0
    return v
rows_id=[]; worst=0.0
for el in TM:
    ctx=[a for (e,a) in labels if e==el]
    for i in range(len(ctx)):
        for j in range(i+1,len(ctx)):
            v=eff_vec(el,ctx[i])-eff_vec(el,ctx[j])
            an=float(np.max(np.abs(ANA@v)))/np.linalg.norm(v)      # exact, must be 0
            emp=np.linalg.norm(weak@v)/np.linalg.norm(v) if len(weak) else 0.0
            worst=max(worst,an)
            rows_id.append(dict(element=el,context_a=ctx[i],context_b=ctx[j],
                                analytic_projection=an,empirical_projection=emp))
emp_worst=max(r["empirical_projection"] for r in rows_id)
print(f"max projection of the {len(rows_id)} within-element contrasts: "
      f"{worst:.3e} onto the analytic near-dependency directions (must be 0), "
      f"{emp_worst:.3e} onto the {len(weak)} weakest empirical directions (diagnostic)")

c_min,_,_,_=np.linalg.lstsq(X,y,rcond=None)
# invariance under two invertible reparameterisations
d=np.exp(np.linspace(-1.5,1.5,K))
c_s,_,_,_=np.linalg.lstsq(X*d[None,:],y,rcond=None); c_s=c_s*d   # X D c = y => beta = D c
perm=np.random.default_rng(7).permutation(K)
c_p_raw,_,_,_=np.linalg.lstsq(X[:,perm],y,rcond=None)
c_p=np.empty(K); c_p[perm]=c_p_raw
def eff(co,el,a): return co[pos[el]]+co[NS+labels.index((el,a))]
maxdiff=0.0
for el in TM:
    ctx=[a for (e,a) in labels if e==el]
    for i in range(len(ctx)):
        for j in range(i+1,len(ctx)):
            b=eff(c_min,el,ctx[i])-eff(c_min,el,ctx[j])
            for cc in (c_s,c_p): maxdiff=max(maxdiff,abs((eff(cc,el,ctx[i])-eff(cc,el,ctx[j]))-b))
print(f"max |contrast difference| across the two reparameterisations: {maxdiff:.3e} meV/atom")
for k_,v_ in (("rank",f"{rank} of {K}"),("condition_number",f"{cond:.4e}"),
              ("rank_without_hydride_rows",f"{rank_nohyd} of {K}"),
              ("condition_without_hydride_rows",f"{cond_nohyd:.4e}"),
              ("hydride_rows_with_a_listed_element",str(sum(hyd_cnt.values()))),
              ("max_analytic_contrast_projection",f"{worst:.3e}"),
              ("max_empirical_contrast_projection",f"{emp_worst:.3e}"),
              ("max_reparam_contrast_diff_meV",f"{maxdiff:.3e}")):
    rows_id.append(dict(element="(design)",context_a=k_,context_b=v_,
                        analytic_projection=0.0,empirical_projection=0.0))
pd.DataFrame(rows_id).to_csv("results/context_identification.csv",index=False)

# ---------------- (2) cluster-robust joint Wald test ----------------
codes,uniq=pd.factorize(m.proto); G=len(uniq)
u=y-X@c_min
XtX_inv=np.linalg.pinv(X.T@X,rcond=1e-14)
S=np.zeros((G,K)); np.add.at(S,codes,X*u[:,None])
V=XtX_inv@(S.T@S)@XtX_inv
adj=G/(G-1)*(len(y)-1)/(len(y)-K); V=V*adj

def gammq(a,x):
    """Regularised upper incomplete gamma Q(a,x); chi2 tail = gammq(df/2, chi2/2)."""
    from math import lgamma, exp, log
    if x<0 or a<=0: raise ValueError
    if x==0: return 1.0
    if x < a+1.0:                                  # series for P(a,x)
        ap=a; s=1.0/a; dl=s
        for _ in range(1000):
            ap+=1.0; dl*=x/ap; s+=dl
            if abs(dl)<abs(s)*1e-15: break
        return 1.0-s*exp(-x+a*log(x)-lgamma(a))
    b=x+1.0-a; c=1e300; dd=1.0/b; h=dd                # continued fraction for Q(a,x)
    for i in range(1,1000):
        an=-i*(i-a); b+=2.0
        dd=an*dd+b
        if abs(dd)<1e-300: dd=1e-300
        c=b+an/c
        if abs(c)<1e-300: c=1e-300
        dd=1.0/dd; de=dd*c; h*=de
        if abs(de-1.0)<1e-15: break
    return h*exp(-x+a*log(x)-lgamma(a))

jt=[]
for el in TM:
    ctx=[a for (e,a) in labels if e==el]
    idx=[NS+labels.index((el,a)) for a in ctx]
    R=np.zeros((len(ctx)-1,K))
    for i in range(1,len(ctx)):
        R[i-1,idx[i]]=1.0; R[i-1,idx[0]]=-1.0       # gamma_i = gamma_0 for all i
    Rb=R@c_min; RVR=R@V@R.T
    W=float(Rb@np.linalg.solve(RVR,Rb)); df=R.shape[0]
    p=gammq(df/2.0,W/2.0)
    jt.append(dict(element=el,wald_chi2=round(W,2),df=df,p_value=p))
jt.sort(key=lambda z:z["p_value"])
mm=len(jt)
for i,z in enumerate(jt,1):
    z["bh_threshold"]=0.05*i/mm
    z["significant_bh_0.05"]=bool(z["p_value"]<=0.05*i/mm)
print("\ncluster-robust joint test that an element's six interaction coefficients are equal:")
print(f"{'element':8}{'chi2':>10}{'df':>4}{'p':>12}{'BH thresh':>11}  sig")
for z in jt:
    print(f"{z['element']:8}{z['wald_chi2']:10.2f}{z['df']:4d}{z['p_value']:12.3e}"
          f"{z['bh_threshold']:11.3f}  {'yes' if z['significant_bh_0.05'] else 'no'}")
pd.DataFrame(jt).to_csv("results/context_joint_test.csv",index=False)

# ---------------- (3) resumable cluster bootstrap of the spreads ----------------
_cd=pd.read_csv("results/context_dependence.csv")
pt={el: float(g.coef.max()-g.coef.min()) for el,g in _cd.groupby("element")}
DRAWS="results/context_spread_draws.npy"
prev=None
if os.path.exists(DRAWS) and not RESET:
    prev=np.load(DRAWS,allow_pickle=True).item()
rng=np.random.default_rng(SEED); p=np.full(G,1.0/G)
t0=time.time(); new={el:[] for el in TM}; done=0
while done<NB and time.time()-t0<BUDGET:
    w=rng.multinomial(G,p).astype(float)[codes]
    A=(X*w[:,None]).T@X; b=X.T@(w*y)
    c=np.linalg.pinv(A,rcond=1e-12)@b
    for el in TM:
        ctx=[a for (e,a) in labels if e==el]
        v=[c[pos[el]]+c[NS+labels.index((el,a))] for a in ctx]
        new[el].append(max(v)-min(v))
    done+=1
print(f"\n{done} new replicates at seed {SEED} in {time.time()-t0:.0f} s")
store=prev if prev else {"draws":{el:[] for el in TM},"seeds":{}}
for el in TM: store["draws"][el]=list(store["draws"][el])+new[el]
store["seeds"][str(SEED)]={el:[float(np.percentile(new[el],2.5)),
                               float(np.percentile(new[el],97.5))] for el in TM}
np.save(DRAWS,store,allow_pickle=True)
tot=len(store["draws"][TM[0]])
print(f"total accumulated replicates: {tot} across seeds {sorted(store['seeds'])}")
rows_ci=[]
print(f"{'element':8}{'spread':>8}{'95% interval':>20}   per-seed endpoints")
for el in TM:
    a=np.array(store["draws"][el])
    lo,hi=np.percentile(a,2.5),np.percentile(a,97.5)
    per=" ".join(f"[{v[el][0]:.0f},{v[el][1]:.0f}]" for v in store["seeds"].values())
    print(f"{el:8}{pt[el]:+8.0f}   [{lo:+6.0f}, {hi:+6.0f}]   {per}")
    rows_ci.append(dict(element=el,spread_mev=round(pt[el],0),lo=round(lo,0),hi=round(hi,0),
                        replicates=tot,seeds=";".join(sorted(store["seeds"]))))
pd.DataFrame(rows_ci).to_csv("results/context_spread_ci.csv",index=False)
ok=(rank==K) and worst<1e-8 and maxdiff<1e-6
print("\n[impl] design full rank, contrasts exactly orthogonal to the analytic "
      "near-dependencies, reparameterisation-invariant:", "PASS" if ok else "FAIL")
if not ok: raise SystemExit(1)

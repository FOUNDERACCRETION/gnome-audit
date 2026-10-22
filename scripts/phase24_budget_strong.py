"""The harder version of the budget test.

The paper currently compares composition against hull distance used as a raw ranking.
A fair competitor is hull distance MODELLED: if a flexible function of the PBE stability
margin already captures what composition is proxying for, the combined gain should
collapse. Policies, all evaluated on identical held-out clusters:

  random        analytic expectation b
  composition   cell rate fitted on the estimation half
  hull_rank     closest to the hull first, no fitting (the paper's current baseline)
  hull_model    logistic on 20 quantile-bin dummies of hull distance -- effectively a
                nonparametric fit in hull distance alone
  combined      hull bins + cell rate
  combined_cf   as combined, but cell rates cross-fitted 5-fold by cluster within the
                estimation half, so the rate entering the model is never fitted on the
                material it scores

Also: do a few small high-rate cells drive the gain? Re-evaluate with the top-3 cells by
fitted rate removed from the ranking advantage.
"""
import ast, re
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
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
def irls(X,y,it=25):
    b=np.zeros(X.shape[1])
    for _ in range(it):
        p=1/(1+np.exp(-X@b)); W=np.clip(p*(1-p),1e-9,None)
        try: b=b+np.linalg.solve((X*W[:,None]).T@X+1e-6*np.eye(X.shape[1]),X.T@(y-p))
        except np.linalg.LinAlgError: break
    return b
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C]).rename(columns={C:"pbe"})
s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(int)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
m=r.merge(s,on="MaterialId",how="inner"); m=m[m.pbe<=0].copy(); m["hull"]=-m.pbe
key=["n_el","anion","Ln","d3"]
protos=m.proto.unique()
print(f"eligible {len(m)}, clusters {len(protos)}")
# fixed quantile knots for the hull basis, from the full eligible set
edges=np.quantile(m.hull.values,np.linspace(0,1,21)); edges[0]-=1e-9; edges[-1]+=1e-9
def hull_basis(v):
    b=np.digitize(v,edges[1:-1])
    X=np.zeros((len(v),20))
    X[np.arange(len(v)),b]=1.0
    return X[:,1:]            # drop first bin as reference
BUD=(0.05,0.10,0.20,0.30,0.50); REP=20
POL=("random","composition","hull_rank","hull_model","combined","combined_cf")
acc={k:{b:[] for b in BUD} for k in POL}
for rep in range(REP):
    perm=RNG.permutation(protos); cut=len(perm)//2
    fitset=set(perm[:cut]); mask=m.proto.isin(fitset)
    fit,hold=m[mask],m[~mask].copy()
    if hold.flip.sum()==0: continue
    rates=fit.groupby(key).flip.mean().rename("rate").reset_index()
    g=fit.flip.mean()
    hold=hold.merge(rates,on=key,how="left"); hold["rate"]=hold["rate"].fillna(g)
    f2=fit.merge(rates,on=key,how="left"); f2["rate"]=f2["rate"].fillna(g)
    # cross-fitted cell rates inside the estimation half
    fp=fit.proto.unique(); fold=dict(zip(fp,RNG.integers(0,5,len(fp))))
    fit=fit.assign(_f=fit.proto.map(fold))
    cf=np.empty(len(fit))
    for k in range(5):
        tr=fit[fit._f!=k]; rt=tr.groupby(key).flip.mean()
        te=fit._f.values==k
        cf[te]=[rt.get(tuple(v),tr.flip.mean()) for v in fit.loc[te,key].itertuples(index=False)]
    Hf=hull_basis(fit.hull.values); Hh=hull_basis(hold.hull.values)
    one_f=np.ones((len(fit),1)); one_h=np.ones((len(hold),1))
    yf=fit.flip.values.astype(float)
    b_h=irls(np.hstack([one_f,Hf]),yf)
    b_c=irls(np.hstack([one_f,Hf,f2.rate.values[:,None]]),yf)
    b_x=irls(np.hstack([one_f,Hf,cf[:,None]]),yf)
    sc={"composition":hold.rate.values,
        "hull_rank":-hold.hull.values,
        "hull_model":(np.hstack([one_h,Hh])@b_h),
        "combined":(np.hstack([one_h,Hh,hold.rate.values[:,None]])@b_c),
        "combined_cf":(np.hstack([one_h,Hh,hold.rate.values[:,None]])@b_x)}
    tot=hold.flip.sum(); n=len(hold); tie=RNG.permutation(n)
    for name,v in sc.items():
        order=np.lexsort((tie,-v)); fl=hold.flip.values[order]
        for bgt in BUD: acc[name][bgt].append(fl[:int(n*bgt)].sum()/tot)
    for bgt in BUD: acc["random"][bgt].append(int(n*bgt)/n)
print(f"\nFlips caught, mean over {REP} cluster-wise splits")
print(f"{'budget':>7}"+"".join(f"{p:>14}" for p in POL))
out=[]
for bgt in BUD:
    rec={"budget":bgt}; line=f"{bgt*100:6.0f}%"
    for p in POL:
        a=np.array(acc[p][bgt]); rec[p]=a.mean(); rec[p+"_sd"]=a.std(ddof=1) if len(a)>1 else 0
        line+=f"{a.mean()*100:13.1f}%"
    out.append(rec); print(line)
pd.DataFrame(out).to_csv("results/budget_strong.csv",index=False)
print("\nPaired per-split differences at a 20% budget (the question that matters):")
_pr=[]
for a,b in (("hull_model","hull_rank"),("combined","hull_model"),("combined_cf","hull_model")):
    d=np.array(acc[a][0.20])-np.array(acc[b][0.20])
    print(f"  {a:12} - {b:12} {d.mean()*100:+6.2f} pp  (sd {d.std(ddof=1)*100:.2f}, "
          f"wins {int((d>0).sum())}/{len(d)})")
    _pr.append(dict(contrast=a+"-"+b,mean_pp=round(d.mean()*100,3),
                    sd_pp=round(d.std(ddof=1)*100,3),wins=int((d>0).sum()),n=len(d)))
pd.DataFrame(_pr).to_csv("results/budget_strong_paired.csv",index=False)
print("wrote results/budget_strong_paired.csv")

"""Referee round 4, item 2: the budget model benchmarked chemistry-ranking only against
uniform random. The obvious competitor is the PBE decomposition energy itself: a
maintainer already has it for every candidate, it is free, and Section 4.4 shows flips
concentrate near the hull. If hull distance alone wins, the paper's operational
recommendation is wrong.

Ranked policies compared on held-out materials, cluster-wise splits so that structural
prototypes never straddle fit and evaluation:
  random          uniform sampling
  chemistry       cell flip rate fitted on the fit half (the paper's current policy)
  hull            closest to the PBE hull first; no chemistry, no fitting
  combined        logistic on cell rate + hull distance, fitted on the fit half
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

s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C])
s=s.rename(columns={C:"pbe"}); s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(int)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
r=r.merge(s,on="MaterialId",how="inner")
r=r[r.pbe<=0].copy()                      # must be PBE-stable to flip
r["hull"]=-r.pbe                          # distance below the hull, positive
print(f"materials with a PBE energy and PBE-stable: {len(r)}  flip {r.flip.mean()*100:.2f}%")
key=["n_el","anion","Ln","d3"]
protos=r.proto.unique()

def logistic_fit(X,y,iters=25):
    X=np.column_stack([np.ones(len(X)),X]); b=np.zeros(X.shape[1])
    for _ in range(iters):
        p=1/(1+np.exp(-X@b)); W=np.clip(p*(1-p),1e-9,None)
        try: b=b+np.linalg.solve((X*W[:,None]).T@X, X.T@(y-p))
        except np.linalg.LinAlgError: break
    return b

BUD=(0.05,0.10,0.20,0.30,0.50); REP=20
acc={k:{b:[] for b in BUD} for k in ("random","chemistry","hull","combined")}
for rep in range(REP):
    perm=RNG.permutation(protos); cut=len(perm)//2
    fitset=set(perm[:cut])
    m=r.proto.isin(fitset)
    fit, hold = r[m], r[~m].copy()
    if hold.flip.sum()==0: continue
    rates=fit.groupby(key).flip.mean().rename("rate").reset_index()
    hold=hold.merge(rates,on=key,how="left")
    hold["rate"]=hold["rate"].fillna(fit.flip.mean())
    f2=fit.merge(rates,on=key,how="left"); f2["rate"]=f2["rate"].fillna(fit.flip.mean())
    b=logistic_fit(np.column_stack([f2.rate.values,f2.hull.values]),f2.flip.values)
    hold["comb"]=1/(1+np.exp(-(b[0]+b[1]*hold.rate.values+b[2]*hold.hull.values)))
    tot=hold.flip.sum(); n=len(hold)
    # Break ties with a seeded random key rather than input order: where the budget
    # cuts through a tied cell, which members are taken changes the realised catch.
    jitter=RNG.permutation(len(hold))
    hold=hold.assign(_tie=jitter)
    ords={"chemistry":hold.sort_values(["rate","_tie"],ascending=[False,True]).flip.values,
          "hull":hold.sort_values(["hull","_tie"],ascending=[True,True]).flip.values,
          "combined":hold.sort_values(["comb","_tie"],ascending=[False,True]).flip.values}
    for bgt in BUD:
        k=int(n*bgt)
        for name,v in ords.items(): acc[name][bgt].append(v[:k].sum()/tot)
        acc["random"][bgt].append(k/n)

# paired per-split differences: the policies share each split, so the marginal s.d.s
# are not the uncertainty on a difference between them.
print("\nPaired per-split differences (the correct uncertainty for policy comparisons):")
for a,b in (("hull","chemistry"),("combined","hull"),("combined","chemistry")):
    for bgt in BUD:
        d=np.array(acc[a][bgt])-np.array(acc[b][bgt])
        if bgt==0.20:
            print(f"  {a} - {b} at 20% budget: {d.mean()*100:+.2f} pp "
                  f"(paired s.d. {d.std(ddof=1)*100:.2f}, min {d.min()*100:+.2f}, "
                  f"max {d.max()*100:+.2f}, wins {int((d>0).sum())}/{len(d)})")
print()
print(f"Flips caught, mean over {REP} cluster-wise splits (s.d. in brackets)")
print(f"{'budget':>7}{'random':>14}{'chemistry':>18}{'hull dist':>18}{'combined':>18}")
rows=[]
for bgt in BUD:
    line=f"{bgt*100:6.0f}%"
    rec={"budget":bgt}
    for name in ("random","chemistry","hull","combined"):
        a=np.array(acc[name][bgt]); line+=f"{a.mean()*100:11.1f} [{a.std(ddof=1)*100:.1f}]"
        rec[name]=a.mean(); rec[name+"_sd"]=a.std(ddof=1)
    rows.append(rec); print(line)
out=pd.DataFrame(rows); out.to_csv("results/budget_baselines.csv",index=False)
print("\n-> results/budget_baselines.csv")
d=out[out.budget==0.20].iloc[0]
print(f"\nAT A 20% BUDGET: chemistry {d.chemistry*100:.1f}%, hull distance {d.hull*100:.1f}%, "
      f"combined {d.combined*100:.1f}%, random {d.random*100:.1f}%")
best=max(("chemistry",d.chemistry),("hull",d.hull),("combined",d.combined),key=lambda x:x[1])
print(f"BEST SINGLE POLICY AT 20%: {best[0]}")

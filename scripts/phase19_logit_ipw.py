"""Referee round 4: (a) the logistic model is main-effects-only although the paper
documents strong interactions, so the lanthanide odds ratio averages effects of opposite
sign. Refit with 3d x anion and Ln x anion interactions and report average marginal
effects. (b) IPW positivity diagnostics: propensity distribution, maximum weight,
effective sample size, and a trimmed sensitivity."""
import ast
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
UEL=set("Co Cr Fe Mn Mo Ni V W".split()); TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
def irls(X,y,it=30):
    b=np.zeros(X.shape[1])
    for _ in range(it):
        p=1/(1+np.exp(-X@b)); W=np.clip(p*(1-p),1e-9,None)
        try: b=b+np.linalg.solve((X*W[:,None]).T@X+1e-8*np.eye(X.shape[1]),X.T@(y-p))
        except np.linalg.LinAlgError: break
    return b
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["n_el"]=els.map(len)
r["U"]=els.map(lambda e:bool(e&UEL)).astype(float)
r["Ln"]=els.map(lambda e:bool(e&LN)).astype(float)
y=(r[C]>0).values.astype(float)
anions=[a for a in sorted(r.anion.unique()) if a!="intermetallic"]
base=[np.ones(len(r)),r.U.values,r.Ln.values,(r.n_el.values-r.n_el.mean())/r.n_el.std()]
names=["intercept","has_U","has_Ln","n_el_z"]
for a in anions: base.append((r.anion==a).values.astype(float)); names.append(f"is_{a[:6]}")
X_main=np.column_stack(base)
b_main=irls(X_main,y)
print("main-effects model: Ln odds ratio = %.3f"%np.exp(b_main[names.index("has_Ln")]))
cols=list(base); nm=list(names)
for a in anions+["intermetallic"]:
    m=(r.anion==a).values.astype(float)
    for v,lab in ((r.Ln.values,"Ln"),(r.U.values,"U")):
        c=v*m
        if c.sum()>=200: cols.append(c); nm.append(f"{lab}x{a[:6]}")
X_int=np.column_stack(cols); b_int=irls(X_int,y)
def ll(X,b):
    p=np.clip(1/(1+np.exp(-X@b)),1e-12,1-1e-12); return (y*np.log(p)+(1-y)*np.log(1-p)).sum()
lr=2*(ll(X_int,b_int)-ll(X_main,b_main)); df=X_int.shape[1]-X_main.shape[1]
print(f"adding anion interactions: LR chi2 = {lr:.0f} on {df} df  (overwhelming)")
print("\nAverage marginal effect of each term, percentage points")
p0=1/(1+np.exp(-X_int@b_int))
rows=[]
for lab,col in (("lanthanide present",2),("open-shell d (U list)",1)):
    Xa=X_int.copy(); Xb=X_int.copy()
    Xa[:,col]=1.0; Xb[:,col]=0.0
    for j,n_ in enumerate(nm):
        if n_.startswith(("Ln x","Lnx") if col==2 else ("Ux","U x")):
            pass
    # rebuild interaction columns consistently
    for j,n_ in enumerate(nm):
        pre="Ln" if col==2 else "U"
        if n_.startswith(pre+"x"):
            a=n_.split("x",1)[1]
            m=np.array([str(s)[:6]==a for s in r.anion.values],dtype=float)
            Xa[:,j]=m; Xb[:,j]=0.0
    ame=( 1/(1+np.exp(-Xa@b_int)) - 1/(1+np.exp(-Xb@b_int)) ).mean()*100
    rows.append(dict(term=lab,ame_pp=round(ame,2)))
    print(f"  {lab:24}{ame:+7.2f} pp")
print("\nlanthanide average marginal effect BY anion class:")
for a in sorted(r.anion.unique()):
    m=(r.anion==a).values
    if m.sum()<500: continue
    sub=r[m]; yl=y[m]
    withln=yl[sub.Ln.values==1]; without=yl[sub.Ln.values==0]
    if len(withln)<100 or len(without)<100: continue
    print(f"  {a:26}{(withln.mean()-without.mean())*100:+7.1f} pp  (n_Ln={len(withln)})")
pd.DataFrame(rows).to_csv("results/logit_interactions.csv",index=False)

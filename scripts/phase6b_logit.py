"""Logistic regression with ANALYTIC standard errors from the Fisher information.
Replaces the under-converged bootstrap in phase6."""
import ast
import numpy as np, pandas as pd
COL="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
UEL=set("Co Cr Fe Mn Mo Ni V W".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["Elements","NSites",COL])
r[COL]=pd.to_numeric(r[COL],errors="coerce"); r=r.dropna(subset=[COL])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
y=(r[COL]>0).values.astype(float)
F=pd.DataFrame({
 "has_U":els.map(lambda e:bool(e&UEL)).astype(float),
 "has_3d_nonU":els.map(lambda e:bool(e&TM3D) and not bool(e&UEL)).astype(float),
 "has_Ln":els.map(lambda e:bool(e&LN)).astype(float),
 "is_oxide":(r.anion=="oxide").astype(float),
 "is_pnictide":(r.anion=="pnictide").astype(float),
 "is_halide":(r.anion=="halide").astype(float),
 "is_chalcogenide":(r.anion=="chalcogenide").astype(float),
 "is_bcs":(r.anion=="boride/carbide/silicide").astype(float),
 "n_elements_z":((r.n_el-r.n_el.mean())/r.n_el.std()).astype(float),
 "log_nsites_z":((np.log(r.NSites)-np.log(r.NSites).mean())/np.log(r.NSites).std()).astype(float),
})
X=np.c_[np.ones(len(F)),F.values]; names=["intercept"]+list(F.columns)
# Newton-Raphson / IRLS to convergence
w=np.zeros(X.shape[1])
for it in range(60):
    z=X@w; p=1/(1+np.exp(-z)); W=p*(1-p)
    g=X.T@(y-p); H=(X*W[:,None]).T@X
    step=np.linalg.solve(H+1e-9*np.eye(len(w)),g); w+=step
    if np.max(np.abs(step))<1e-10: break
print(f"IRLS converged in {it+1} iterations, max|step|={np.max(np.abs(step)):.2e}")
z=X@w; p=1/(1+np.exp(-z)); W=p*(1-p)
cov=np.linalg.inv((X*W[:,None]).T@X); se=np.sqrt(np.diag(cov))
res=pd.DataFrame(dict(feature=names,coef=w,se=se,odds_ratio=np.exp(w),
    or_lo=np.exp(w-1.96*se),or_hi=np.exp(w+1.96*se),z=w/se))
res.to_csv("results/logistic_coefficients.csv",index=False)
print(res.to_string(index=False,float_format=lambda x:f"{x:.4f}"))
ll=np.sum(y*np.log(np.clip(p,1e-12,1))+(1-y)*np.log(np.clip(1-p,1e-12,1)))
pbar=y.mean(); ll0=len(y)*(pbar*np.log(pbar)+(1-pbar)*np.log(1-pbar))
print(f"\nn={len(y)}  McFadden pseudo-R2 = {1-ll/ll0:.4f}")

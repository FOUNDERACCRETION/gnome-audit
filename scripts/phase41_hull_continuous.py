"""Continuous nonlinear adjustment for distance below the PBE hull.

Section 4.9 stratifies by PBE decomposition-energy band, which leaves residual variation in
distance inside each band. This fits the chemistry model with a natural cubic spline in the
continuous PBE decomposition energy instead, on the matched subset where that quantity is
observed, with cluster-robust standard errors over structural prototypes. It also fits a
chemistry-by-depth interaction, so that the question "does the chemistry contrast merely
stand in for hull proximity?" is answered without binning.

Outputs results/hull_continuous.csv and prints an unadjusted/adjusted comparison.
"""
import ast, re
import numpy as np, pandas as pd
CR="Decomposition Energy Per Atom MP"      # r2SCAN decomposition energy (the outcome side)
CP="Decomposition Energy Per Atom MP"       # PBE decomposition energy, SAME reference hull
# as the outcome. The summary file also carries a plain "Decomposition Energy Per Atom"
# column; using it here would define a different eligible population (116,988 PBE-stable
# rather than the 117,018 of Section 3.4) and would mix reference sets across the two sides.
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
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))

print("loading r2SCAN file...",flush=True)
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","NSites","Space Group Number",CR])
r[CR]=pd.to_numeric(r[CR],errors="coerce"); r=r.dropna(subset=[CR])
print("loading summary file for the PBE covariate...",flush=True)
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",CP])
s["pbe_dec"]=pd.to_numeric(s[CP],errors="coerce")
s=s[["MaterialId","pbe_dec"]].dropna(subset=["pbe_dec"])
assert not s.MaterialId.duplicated().any(), "summary MaterialId is not unique"
d=r.merge(s,on="MaterialId",how="inner").reset_index(drop=True)
print(f"r2SCAN rows {len(r)}, matched with an observed PBE decomposition energy {len(d)}")
d=d[d.pbe_dec<=0].reset_index(drop=True)
print(f"restricted to observed PBE-stable (decomposition energy <= 0): {len(d)}")

els=d["Elements"].map(ast.literal_eval).map(set)
d["n_el"]=els.map(len); d["anion"]=els.map(ac)
d["proto"]=d["Space Group Number"].astype(str)+"|"+d.Composition.map(anon)
y=(d[CR]>0).values.astype(float)
dep=d.pbe_dec.values                                  # <= 0; more negative = deeper below hull

def ns_basis(x,df=5):
    """Natural cubic spline basis on x with knots at equally spaced quantiles.
    Truncated-power cubic basis with linear tails imposed beyond the boundary knots."""
    qs=np.linspace(0,1,df+1)
    k=np.quantile(x,qs)
    k=np.unique(k)
    if len(k)<4: raise SystemExit("too few distinct knots")
    kb0,kbK=k[0],k[-1]; inner=k[1:-1]
    def cub(u,t): 
        v=u-t; return np.where(v>0,v**3,0.0)
    # Natural-spline constraint: d(x) = (cub(x,t) - cub(x,kbK))/(kbK-t) differences
    cols=[x]
    dK=(cub(x,inner[-1])-cub(x,kbK))/(kbK-inner[-1]) if len(inner) else None
    for t in inner[:-1]:
        dt=(cub(x,t)-cub(x,kbK))/(kbK-t)
        cols.append(dt-dK)
    B=np.column_stack(cols)
    B=(B-B.mean(0))/np.where(B.std(0)>0,B.std(0),1.0)
    return B,k

def irls(X,y):
    w=np.zeros(X.shape[1])
    for it in range(80):
        p=1/(1+np.exp(-np.clip(X@w,-30,30))); W=np.clip(p*(1-p),1e-10,None)
        g=X.T@(y-p); H=(X*W[:,None]).T@X
        step=np.linalg.solve(H+1e-8*np.eye(len(w)),g); w+=step
        if np.max(np.abs(step))<1e-10: break
    return w,it+1

def cluster_se(X,y,w,codes):
    p=1/(1+np.exp(-np.clip(X@w,-30,30))); W=np.clip(p*(1-p),1e-10,None)
    A=np.linalg.inv((X*W[:,None]).T@X)
    U=X*(y-p)[:,None]
    G=codes.max()+1
    S=np.zeros((G,X.shape[1]))
    np.add.at(S,codes,U)
    meat=S.T@S
    V=A@meat@A
    nG=G; k=X.shape[1]; n=len(y)
    adj=nG/(nG-1)*(n-1)/(n-k)
    return np.sqrt(np.diag(V)*adj)

FEAT=pd.DataFrame({
 "has_U":els.map(lambda e:bool(e&UEL)).astype(float).values,
 "has_3d_nonU":els.map(lambda e:bool(e&TM3D) and not bool(e&UEL)).astype(float).values,
 "has_Ln":els.map(lambda e:bool(e&LN)).astype(float).values,
 "is_oxide":(d.anion=="oxide").astype(float).values,
 "is_pnictide":(d.anion=="pnictide").astype(float).values,
 "is_halide":(d.anion=="halide").astype(float).values,
 "is_chalcogenide":(d.anion=="chalcogenide").astype(float).values,
 "is_bcs":(d.anion=="boride/carbide/silicide").astype(float).values,
 "n_elements_z":((d.n_el-d.n_el.mean())/d.n_el.std()).values,
 "log_nsites_z":((np.log(d.NSites)-np.log(d.NSites).mean())/np.log(d.NSites).std()).values,
})
codes,_=pd.factorize(d.proto); codes=codes.astype(np.int64)
print(f"prototype clusters in this subset: {codes.max()+1}")

rows=[{"model":"population","feature":"n_matched_with_pbe_decomp","odds_ratio":len(d),
        "or_lo":"","or_hi":""},
      {"model":"population","feature":"n_prototype_clusters","odds_ratio":int(codes.max()+1),
        "or_lo":"","or_hi":""}]
def fit(label,extra=None,extra_names=()):
    cols=[np.ones(len(d)),FEAT.values]
    if extra is not None: cols.append(extra)
    X=np.column_stack(cols)
    names=["intercept"]+list(FEAT.columns)+list(extra_names)
    w,nit=irls(X,y); se=cluster_se(X,y,w,codes)
    print(f"\n--- {label} (IRLS {nit} iters, k={X.shape[1]}) ---")
    for nm,b,s_ in zip(names,w,se):
        if nm.startswith(("spline","dep")): continue
        print(f"  {nm:16s} OR {np.exp(b):7.3f}  [{np.exp(b-1.96*s_):6.3f}, {np.exp(b+1.96*s_):6.3f}]")
        rows.append(dict(model=label,feature=nm,odds_ratio=round(float(np.exp(b)),3),
                         or_lo=round(float(np.exp(b-1.96*s_)),3),
                         or_hi=round(float(np.exp(b+1.96*s_)),3)))
    p=1/(1+np.exp(-np.clip(X@w,-30,30)))
    ll=np.sum(y*np.log(np.clip(p,1e-12,1))+(1-y)*np.log(np.clip(1-p,1e-12,1)))
    pb=y.mean(); ll0=len(y)*(pb*np.log(pb)+(1-pb)*np.log(1-pb))
    print(f"  McFadden pseudo-R2 {1-ll/ll0:.4f}")
    rows.append(dict(model=label,feature="mcfadden_r2",odds_ratio=round(float(1-ll/ll0),4),
                     or_lo="",or_hi=""))
    return w,se,X,names

fit("no hull-distance term")
B,knots=ns_basis(dep,df=6)
print(f"\nnatural spline knots (eV/atom): {np.round(knots,4).tolist()}")
w2,se2,X2,n2=fit("natural spline in PBE decomposition energy",
                 B,[f"spline{i+1}" for i in range(B.shape[1])])
# chemistry-by-depth interaction: does the chemistry effect vary with depth?
depz=(dep-dep.mean())/dep.std()
INT=np.column_stack([FEAT["has_U"].values*depz,FEAT["is_oxide"].values*depz])
w3,se3,X3,n3=fit("spline plus has_U and oxide interactions with depth",
                 np.column_stack([B,INT]),
                 [f"spline{i+1}" for i in range(B.shape[1])]+["has_U_x_depth","is_oxide_x_depth"])
for nm in ("has_U_x_depth","is_oxide_x_depth"):
    i=n3.index(nm)
    print(f"  {nm:18s} coef {w3[i]:+.4f}  z={w3[i]/se3[i]:+.2f}")
    rows.append(dict(model="interaction terms",feature=nm,
                     odds_ratio=round(float(np.exp(w3[i])),3),
                     or_lo=round(float(np.exp(w3[i]-1.96*se3[i])),3),
                     or_hi=round(float(np.exp(w3[i]+1.96*se3[i])),3)))
pd.DataFrame(rows).to_csv("results/hull_continuous.csv",index=False)
print("\nwrote results/hull_continuous.csv")

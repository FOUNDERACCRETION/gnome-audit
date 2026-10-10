"""Referee round 2:
(a) cluster-robust sandwich SEs for the logistic model
(b) does the headline survive an ALTERNATIVE clustering definition?"""
import ast,re
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
C="Decomposition Energy Per Atom MP"
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

r=pd.read_csv("data/stable_materials_r2scan.csv",
  usecols=["Composition","Elements","NSites","Space Group Number","Crystal System",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
y=(r[C]>0).values.astype(float)
anonf=r.Composition.map(anon)

# three clustering definitions, coarse -> fine
clu={
 "A_anon_only (coarsest)":      anonf.values,
 "B_spacegroup+anon (primary)": (r["Space Group Number"].astype(str)+"|"+anonf).values,
 "C_sg+anon+elementset":        (r["Space Group Number"].astype(str)+"|"+anonf+"|"+
                                 els.map(lambda e:"".join(sorted(e)))).values,
}
print("=== clustering definitions ===")
for k,v in clu.items():
    s=pd.Series(v).value_counts()
    print(f"{k:30s} clusters {len(s):6d}  mean size {len(v)/len(s):6.2f}  max {s.iloc[0]}")

F=pd.DataFrame({
 "has_U":els.map(lambda e:bool(e&UEL)).astype(float).values,
 "has_3d_nonU":els.map(lambda e:bool(e&TM3D) and not bool(e&UEL)).astype(float).values,
 "has_Ln":els.map(lambda e:bool(e&LN)).astype(float).values,
 "is_oxide":(r.anion=="oxide").astype(float).values,
 "is_pnictide":(r.anion=="pnictide").astype(float).values,
 "is_halide":(r.anion=="halide").astype(float).values,
 "is_chalc":(r.anion=="chalcogenide").astype(float).values,
 "is_bcs":(r.anion=="boride/carbide/silicide").astype(float).values,
 "n_el_z":((r.n_el-r.n_el.mean())/r.n_el.std()).values,
 "lognsites_z":((np.log(r.NSites)-np.log(r.NSites).mean())/np.log(r.NSites).std()).values,
})
X=np.c_[np.ones(len(F)),F.values]; names=["intercept"]+list(F.columns)
w=np.zeros(X.shape[1])
for it in range(60):
    p=1/(1+np.exp(-(X@w))); W=p*(1-p)
    H=(X*W[:,None]).T@X
    step=np.linalg.solve(H+1e-10*np.eye(len(w)),X.T@(y-p)); w+=step
    if np.max(np.abs(step))<1e-11: break
p=1/(1+np.exp(-(X@w))); W=p*(1-p); u=(y-p)[:,None]*X
bread=np.linalg.inv((X*W[:,None]).T@X)
se_naive=np.sqrt(np.diag(bread))

print("\n=== logistic model: naive vs CLUSTER-ROBUST sandwich SEs ===")
res=pd.DataFrame({"feature":names,"coef":w,"OR":np.exp(w),"se_naive":se_naive})
for label,g in clu.items():
    df=pd.DataFrame(u); df["g"]=g
    S=df.groupby("g").sum().values
    meat=S.T@S
    V=bread@meat@bread
    se=np.sqrt(np.diag(V))
    res[f"se_{label[0]}"]=se
    res[f"lo_{label[0]}"]=np.exp(w-1.96*se); res[f"hi_{label[0]}"]=np.exp(w+1.96*se)
show=res[["feature","OR","se_naive","se_A","se_B","se_C","lo_B","hi_B"]]
print(show.to_string(index=False,float_format=lambda x:f"{x:.4f}"))
res.to_csv("results/logit_cluster_robust.csv",index=False)
i=names.index("has_U")
print(f"\nopen-shell d odds ratio {np.exp(w[i]):.3f}")
for L in "ABC":
    print(f"  cluster-robust ({L}) 95% CI [{res[f'lo_{L}'][i]:.3f}, {res[f'hi_{L}'][i]:.3f}]"
          f"   SE inflation {res[f'se_{L}'][i]/se_naive[i]:.2f}x")

# (b) headline cells under each clustering definition
print("\n=== headline quaternary contrast under each clustering definition ===")
r["flip"]=y
hi=r[(r.n_el==4)&(r.anion=="oxide")&(~r.Ln)&(r.d3)]
lo=r[(r.n_el==4)&(r.anion=="intermetallic")&(~r.Ln)&(~r.d3)]
for label,g in clu.items():
    r["_g"]=g
    hi=r[(r.n_el==4)&(r.anion=="oxide")&(~r.Ln)&(r.d3)]
    lo=r[(r.n_el==4)&(r.anion=="intermetallic")&(~r.Ln)&(~r.d3)]
    def cb(sub,nb=500):
        grp=[x.flip.values for _,x in sub.groupby("_g")]; K=len(grp)
        b=np.array([np.concatenate([grp[i] for i in RNG.integers(0,K,K)]).mean()
                    for _ in range(nb)])
        return sub.flip.mean(),np.percentile(b,2.5),np.percentile(b,97.5)
    ph,lh,hh=cb(hi); pl,ll,hl=cb(lo)
    print(f"{label:30s} hi {ph*100:5.1f}% [{lh*100:5.1f},{hh*100:5.1f}]   "
          f"lo {pl*100:4.1f}% [{ll*100:4.1f},{hl*100:4.1f}]   disjoint {hh>0 and lh>hl}")

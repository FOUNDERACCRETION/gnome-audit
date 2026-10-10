"""Extended analyses: per-element, band gap, crystal system, size, hull distance,
and a hand-rolled logistic regression (no sklearn in this environment)."""
import ast, json
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
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
def ci(v,nb=1000):
    n=len(v); p=float(np.mean(v)); d=RNG.binomial(n,p,nb)/n
    return p,float(np.percentile(d,2.5)),float(np.percentile(d,97.5))

r=pd.read_csv("data/stable_materials_r2scan.csv",
    usecols=["MaterialId","Elements","NSites","Volume","Density","Crystal System",
             "Space Group Number",COL])
r[COL]=pd.to_numeric(r[COL],errors="coerce"); r=r.dropna(subset=[COL])
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["U"]=els.map(lambda e:bool(e&UEL))
r["flip"]=(r[COL]>0).astype(int)
r["vol_per_site"]=r.Volume/r.NSites

# bandgap lives in the summary file
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId","Bandgap"])
r=r.merge(s,on="MaterialId",how="left")
print(f"rows {len(r)}, bandgap present for {r.Bandgap.notna().sum()}",flush=True)

out={}

# A. per-element
rows=[]
allel=sorted({x for e in els for x in e})
for el in allel:
    m=els.map(lambda e: el in e).values
    if m.sum()<300: continue
    p,lo,hi=ci(r.flip.values[m])
    rows.append(dict(element=el,n=int(m.sum()),flip=p,lo=lo,hi=hi))
pe=pd.DataFrame(rows).sort_values("flip",ascending=False)
pe.to_csv("results/per_element.csv",index=False)
print(f"\n=== PER-ELEMENT (n>=300), top 15 / bottom 10 of {len(pe)} ===",flush=True)
print(pe.head(15).to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)
print(pe.tail(10).to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)

# B. band gap: metal vs non-metal, stratified by 3d
print("\n=== BAND GAP ===",flush=True)
b=r.dropna(subset=["Bandgap"]).copy()
b["metal"]=b.Bandgap<=0.01
rows=[]
for n_el in (3,4):
    for d3 in (True,False):
        for met in (True,False):
            g=b[(b.n_el==n_el)&(b.d3==d3)&(b.metal==met)]
            if len(g)<300: continue
            p,lo,hi=ci(g.flip.values)
            rows.append(dict(n_el=n_el,has_3d=d3,metal=met,n=len(g),flip=p,lo=lo,hi=hi))
bg=pd.DataFrame(rows); bg.to_csv("results/bandgap_split.csv",index=False)
print(bg.to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)
# gap bins
b["gapbin"]=pd.cut(b.Bandgap,[-.01,.01,.5,1.5,3,20])
gb=b.groupby("gapbin",observed=True).agg(n=("flip","size"),flip=("flip","mean")).reset_index()
gb.to_csv("results/bandgap_bins.csv",index=False)
print(gb.to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)

# C. crystal system
print("\n=== CRYSTAL SYSTEM (n>=500) ===",flush=True)
rows=[]
for cs,g in r.groupby("Crystal System"):
    if len(g)<500: continue
    p,lo,hi=ci(g.flip.values)
    rows.append(dict(crystal_system=cs,n=len(g),flip=p,lo=lo,hi=hi))
cs=pd.DataFrame(rows).sort_values("flip",ascending=False)
cs.to_csv("results/crystal_system.csv",index=False)
print(cs.to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)

# D. hull distance bins x anion (tests Note 2's near-hull claim)
print("\n=== PBE-side proxy: r2SCAN decomp distribution by anion ===",flush=True)
r["dbin"]=pd.cut(r[COL],[-np.inf,-0.1,-0.05,-0.02,-0.005,0,np.inf])
hd=r.pivot_table(index="dbin",columns="anion",values="flip",aggfunc="size",observed=True)
hd.to_csv("results/decomp_distribution.csv")
print(hd.to_string(),flush=True)

# E. size effects
print("\n=== SIZE ===",flush=True)
r["sbin"]=pd.cut(r.NSites,[0,10,16,24,32,100])
sz=r.groupby("sbin",observed=True).agg(n=("flip","size"),flip=("flip","mean")).reset_index()
sz.to_csv("results/size_effect.csv",index=False)
print(sz.to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)

# F. logistic regression, hand-rolled (IRLS-free: plain gradient descent)
print("\n=== LOGISTIC REGRESSION (composition features) ===",flush=True)
feat=pd.DataFrame({
  "has_U":r.U.astype(float),
  "has_3d_nonU":(r.d3&~r.U).astype(float),
  "has_Ln":r.Ln.astype(float),
  "is_oxide":(r.anion=="oxide").astype(float),
  "is_pnictide":(r.anion=="pnictide").astype(float),
  "is_halide":(r.anion=="halide").astype(float),
  "is_chalcogenide":(r.anion=="chalcogenide").astype(float),
  "is_bcs":(r.anion=="boride/carbide/silicide").astype(float),
  "n_elements":(r.n_el-r.n_el.mean())/r.n_el.std(),
  "log_nsites":(np.log(r.NSites)-np.log(r.NSites).mean())/np.log(r.NSites).std(),
})
X=np.c_[np.ones(len(feat)),feat.values]; y=r.flip.values.astype(float)
w=np.zeros(X.shape[1]); lr=0.5
for i in range(4000):
    z=X@w; p=1/(1+np.exp(-z)); g=X.T@(p-y)/len(y); w-=lr*g
names=["intercept"]+list(feat.columns)
# bootstrap SEs on coefficients
B=[]
for _ in range(20):
    idx=RNG.integers(0,len(y),len(y)); Xb,yb=X[idx],y[idx]; wb=np.zeros(X.shape[1])
    for i in range(600):
        p=1/(1+np.exp(-(Xb@wb))); wb-=lr*(Xb.T@(p-yb)/len(yb))
    B.append(wb)
B=np.array(B)
lrres=pd.DataFrame(dict(feature=names,coef=w,odds_ratio=np.exp(w),
                        se=B.std(0),lo=np.exp(np.percentile(B,2.5,axis=0)),
                        hi=np.exp(np.percentile(B,97.5,axis=0))))
lrres.to_csv("results/logistic_coefficients.csv",index=False)
print(lrres.to_string(index=False,float_format=lambda x:f"{x:.4f}"),flush=True)
print("\ndone",flush=True)

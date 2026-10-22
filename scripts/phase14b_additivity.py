"""Is Delta additive in composition (bookkeeping) or context-dependent (chemistry)?

Why the naive test fails: the composition cells are defined by element PRESENCE, so an
additive per-element model spans almost the same space as the cell indicators. Removing
the additive part necessarily removes the cell gradient; that is collinearity, not
evidence.

The discriminating question is whether a given element's contribution to Delta depends on
what it is combined with. An additive model says Fe contributes the same per unit atomic
fraction in an oxide and in a boride. Fit that model, then let Fe/Co/Ni/Mn/Cr interact
with anion class, and measure how much the fit improves and how far the per-context
coefficients spread.
"""
import ast, re
import numpy as np, pandas as pd
RNG=np.random.default_rng(20261009)
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
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C]).rename(columns={C:"pbe"})
s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Composition","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).rename(columns={C:"scan"})
els=r["Elements"].map(ast.literal_eval).map(set); r["anion"]=els.map(ac)
m=r.merge(s,on="MaterialId",how="inner"); m=m[m.pbe<=0].copy()
m["delta"]=(m.scan-m.pbe)*1000.0
def frac(c):
    d={}
    for sym,n in re.findall(r"([A-Z][a-z]?)(\d+)",c): d[sym]=d.get(sym,0)+int(n)
    t=sum(d.values()); return {k:v/t for k,v in d.items()}
F=m.Composition.map(frac); syms=sorted({k for d in F for k in d}); pos={s_:i for i,s_ in enumerate(syms)}
X0=np.zeros((len(m),len(syms)))
for i,d in enumerate(F.values):
    for k,v in d.items(): X0[i,pos[k]]=v
y=m.delta.values
def r2(X,y):
    c,_,_,_=np.linalg.lstsq(X,y,rcond=None); p=X@c
    return 1-((y-p)**2).sum()/((y-y.mean())**2).sum(), c
R2a,_=r2(X0,y)
anions=sorted(m.anion.unique()); av=m.anion.values
cols=[]; labels=[]
for el in TM:
    xe=X0[:,pos[el]]
    for a in anions:
        v=xe*(av==a)
        if (v>0).sum()>=100: cols.append(v); labels.append((el,a))
X1=np.column_stack([X0]+cols)
R2b,cb=r2(X1,y)
print(f"additive per-element model      R^2 = {R2a:.4f}")
print(f"  + element x anion interaction R^2 = {R2b:.4f}   (improvement {R2b-R2a:+.4f})")
print(f"  interaction terms: {len(cols)}, n = {len(y)}")
base={el:cb[pos[el]] for el in TM}
print("\nEffective coefficient per unit atomic fraction, by context (meV/atom):")
print(f"{'element':9}{'anion':26}{'coef':>9}")
rows=[]
for (el,a),c in zip(labels,cb[len(syms):]):
    eff=base[el]+c
    rows.append(dict(element=el,anion=a,coef=round(eff,0)))
for el in TM:
    sub=[x for x in rows if x["element"]==el]
    if not sub: continue
    sub.sort(key=lambda x:-x["coef"])
    for x in sub: print(f"{x['element']:9}{x['anion']:26}{x['coef']:+9.0f}")
    sp=sub[0]["coef"]-sub[-1]["coef"]
    print(f"{'':9}{'spread across contexts':26}{sp:+9.0f}\n")
pd.DataFrame(rows).to_csv("results/context_dependence.csv",index=False)
print("An additive model requires every spread above to be zero.")

"""Register items 01/04/05, 30, 32, 37/38, 40, 42, 43."""
import ast, re, json
import numpy as np, pandas as pd
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

r=pd.read_csv("data/stable_materials_r2scan.csv")
GAP=[c for c in r.columns if "gap" in c.lower()]
print("gap-like columns:", GAP)
r[C]=pd.to_numeric(r[C],errors="coerce")
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(int)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)

s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C])
s=s.rename(columns={C:"pbe"}); s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
m=r.merge(s[["MaterialId","pbe"]],on="MaterialId",how="left")
matched=m.pbe.notna(); bad=matched & (m.pbe>0)
print("\n=== A. population reconciliation ===")
print("r2scan rows                          %d"%len(r))
print("matched to summary                   %d"%matched.sum())
print("matched and PBE-stable               %d"%(matched&(m.pbe<=0)).sum())
print("matched and NOT PBE-stable (the 25)  %d"%bad.sum())
cnt=m.groupby(["n_el","anion","Ln","d3"]).size()
prim=set(cnt[cnt>=500].index)
inprim=pd.Series([ (a,b2,c2,d) in prim for a,b2,c2,d in zip(m.n_el,m.anion,m.Ln,m.d3)], index=m.index)
elig_prim = inprim & m.n_el.isin([3,4]) & (m.anion!="hydride")
print("primary-cell population              %d"%elig_prim.sum())
print("  of the 25, in primary cells:       %d"%(bad&elig_prim).sum())
if GAP:
    gc=GAP[0]; g=pd.to_numeric(r[gc],errors="coerce")
    print("gap column: %s"%gc)
    print("  finite positive %d   exactly zero %d   missing %d"%((g>0).sum(),(g==0).sum(),g.isna().sum()))
    print("  of the 25, finite positive gap:  %d"%(bad&(g>0)).sum())
    print("\n=== B. exact-zero gap provenance (item 30) ===")
    z=r.loc[g==0,gc].astype(str)
    print("  raw string forms:", dict(z.value_counts()))
    print("  gaps in (0,1e-3): %d   min positive gap %.6g   negatives %d"%(((g>0)&(g<1e-3)).sum(), g[g>0].min(), (g<0).sum()))
    zc=r.loc[g==0]
    print("  anion mix of the zero-gap rows:", dict(zc.anion.value_counts()))
    print("  flip rate among zero-gap rows: %.2f%%"%(zc.flip.mean()*100))
print("\nfull-file positive-r2SCAN rate          %.4f%%  (n=%d)"%(r.flip.mean()*100,len(r)))
keep=~bad
print("dropping the 25 observed-unstable       %.4f%%  (n=%d)"%(m.loc[keep,'flip'].mean()*100,keep.sum()))
print("primary-cell rate, all                  %.4f%%  (n=%d)"%(m.loc[elig_prim,'flip'].mean()*100,elig_prim.sum()))
print("primary-cell rate, dropping the 25      %.4f%%  (n=%d)"%(m.loc[elig_prim&keep,'flip'].mean()*100,(elig_prim&keep).sum()))

print("\n=== F. budget cluster pool (item 43) ===")
elig=m[matched&(m.pbe<=0)]
print("eligible for budget evaluation: %d"%len(elig))
print("clusters among ALL r2scan rows: %d"%r.proto.nunique())
print("clusters among eligible rows  : %d   half = %d"%(elig.proto.nunique(), elig.proto.nunique()//2))

pd.DataFrame([
  dict(quantity="r2scan_rows",value=len(r)),
  dict(quantity="matched",value=int(matched.sum())),
  dict(quantity="matched_pbe_stable",value=int((matched&(m.pbe<=0)).sum())),
  dict(quantity="matched_not_pbe_stable",value=int(bad.sum())),
  dict(quantity="of_those_in_primary_cells",value=int((bad&elig_prim).sum())),
  dict(quantity="primary_cell_n",value=int(elig_prim.sum())),
  dict(quantity="primary_cell_n_excl_25",value=int((elig_prim&keep).sum())),
  dict(quantity="budget_cluster_pool",value=elig.proto.nunique()),
  dict(quantity="budget_clusters_per_split",value=elig.proto.nunique()//2),
  dict(quantity="all_clusters",value=r.proto.nunique()),
]).to_csv("results/populations.csv",index=False)
print("wrote results/populations.csv")

print("\n=== E. pnictide threshold values (item 42) ===")
tc=pd.read_csv("results/threshold_curves.csv")
cols=["threshold_eV","overall","oxide","pnictide","intermetallic"]
print(tc[cols].to_string(index=False))

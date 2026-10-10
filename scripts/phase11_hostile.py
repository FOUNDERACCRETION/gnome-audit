"""Hostile referee pass. Questions the paper does not currently answer:
 A. Does the chemistry gradient survive WITHIN PBE hull-distance bins? (pre-reg gate 6)
 B. Is r2SCAN shifting certain chemistries systematically, or just adding noise?
 C. Does the ordering hold against the 'All' reference column (25.8% flip)?
 D. Is Fe's 42.8% Fe, or is it Fe-oxides?
 E. Were all recomputed materials actually PBE-stable?"""
import ast
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"; CA="Decomposition Energy Per Atom All"
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
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C])
s=s.rename(columns={C:"pbe"}); s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C,CA])
r[C]=pd.to_numeric(r[C],errors="coerce"); r[CA]=pd.to_numeric(r[CA],errors="coerce")
r=r.rename(columns={C:"scan",CA:"scan_all"})
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r.scan>0).astype(int); r["flip_all"]=(r.scan_all>0).astype(int)
m=r.merge(s,on="MaterialId",how="inner")
print(f"matched with PBE energy: {len(m)} of {len(r)}")

# E. were they PBE-stable?
print(f"\nE. PBE decomp > 0 among recomputed (i.e. NOT PBE-stable): {(m.pbe>0).sum()}"
      f"  ({(m.pbe>0).mean()*100:.2f}%)")
print(f"   PBE decomp distribution: median {m.pbe.median():.4f}, "
      f"10th pct {m.pbe.quantile(.1):.4f}, 90th pct {m.pbe.quantile(.9):.4f}")
m=m[m.pbe<=0].copy()
print(f"   restricting to PBE-stable: {len(m)}")

# A. chemistry gradient within PBE hull-distance bins
m["bin"]=pd.cut(m.pbe,[-np.inf,-0.1,-0.05,-0.02,-0.01,0.0],
                labels=["<-0.1","-0.1..-0.05","-0.05..-0.02","-0.02..-0.01","-0.01..0"])
q=m[m.n_el==4]
print("\nA. QUATERNARY flip rate by PBE hull distance (rows) x anion/3d (cols)")
q["cell"]=q.anion.str[:6]+np.where(q.d3,"+3d","")
pv=q.pivot_table(index="bin",columns="cell",values="flip",aggfunc="mean",observed=True)
cnt=q.pivot_table(index="bin",columns="cell",values="flip",aggfunc="size",observed=True)
pv=pv.where(cnt>=200)
keep=["oxide+3d","oxide","pnicti+3d","pnicti","interm+3d","interm"]
keep=[k for k in keep if k in pv.columns]
print((pv[keep]*100).round(1).to_string())
print("counts:"); print(cnt[keep].to_string())
# within-bin spread
print("\n   within-bin spread (oxide+3d minus intermetallic), pp:")
for b in pv.index:
    if "oxide+3d" in pv.columns and "interm" in pv.columns:
        a,c=pv.loc[b,"oxide+3d"],pv.loc[b,"interm"]
        if pd.notna(a) and pd.notna(c): print(f"     {b:14s} {(a-c)*100:+6.1f}")
# PBE hull-distance distribution by cell: are oxides closer to hull?
print("\n   median PBE decomp energy by cell (are oxides just nearer the hull?):")
print((q.groupby("cell").pbe.median()[keep]*1000).round(1).rename("meV/atom").to_string())

# B. the SHIFT: delta = scan - pbe, by cell
m["delta"]=m.scan-m.pbe
print("\nB. r2SCAN minus PBE decomposition energy (meV/atom), quaternary, by cell")
g=q.assign(delta=q.scan-q.pbe).groupby("cell").delta
tab=pd.DataFrame({"n":g.size(),"mean":g.mean()*1000,"median":g.median()*1000,
                  "sd":g.std()*1000,"frac_up":g.apply(lambda v:(v>0).mean())})
print(tab.loc[keep].round(1).to_string())
print("   (mean>0 = r2SCAN systematically DEstabilises; frac_up = share shifted upward)")

# C. ordering under the 'All' reference column
print("\nC. anion-class flip rate: MP reference vs 'All' reference")
t=r.groupby("anion").agg(n=("flip","size"),flip_MP=("flip","mean"),flip_All=("flip_all","mean"))
t=t[t.n>=500].sort_values("flip_MP",ascending=False)
print((t.assign(flip_MP=t.flip_MP*100,flip_All=t.flip_All*100)).round(1).to_string())
# both series are already ranks, so Pearson on ranks IS Spearman; avoids a scipy dependency
print(f"   rank correlation of the two orderings: "
      f"{t.flip_MP.rank().corr(t.flip_All.rank()):.3f}")

# D. open-shell d element excess WITHIN anion class.
# Per-element rates are confounded with the anion classes an element tends to occur in,
# so the quantity to read is the excess over the no-that-element baseline inside each class.
# This is the source table for Section 4.3 (tab:feexcess).
print("\nD. Within-anion excess flip rate for Fe, Co, Ni")
rows=[]
for sym in ("Fe","Co","Ni"):
    has=els.map(lambda e,s=sym: s in e).values
    with_=r[has]; without=r[~has]
    a=with_.groupby("anion").agg(n=("flip","size"),flip=("flip","mean"))
    b=without.groupby("anion").agg(flip=("flip","mean"))
    j=a.join(b,rsuffix="_without")
    j=j[j.n>=100]
    for anion,v in j.iterrows():
        rows.append(dict(element=sym,anion=anion,n_with=int(v.n),
                         flip_with=round(v.flip*100,1),
                         flip_without=round(v.flip_without*100,1),
                         excess_pp=round((v.flip-v.flip_without)*100,1)))
exc=pd.DataFrame(rows).sort_values(["element","excess_pp"],ascending=[True,False])
print(exc.to_string(index=False))
exc.to_csv("results/within_anion_excess.csv",index=False)
print("   -> results/within_anion_excess.csv")

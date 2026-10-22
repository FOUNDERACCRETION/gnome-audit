"""Referee round 4: (a) the band-gap control collapsed a non-monotonic relationship into a
binary split; redo it in bins crossed with 3d and anion class. (b) the directional-shift
table reported 'share shifted up' rounded to 10% and used a cell definition that differs
by a few materials from the primary cells; recompute on the primary definition with
quartiles and 1% resolution."""
import ast
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
s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C,"Bandgap"]).rename(columns={C:"pbe"})
s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
s["Bandgap"]=pd.to_numeric(s.Bandgap,errors="coerce").replace([np.inf,-np.inf],np.nan)
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).rename(columns={C:"scan"})
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r.scan>0).astype(int)
r=r.merge(s[["MaterialId","Bandgap"]],on="MaterialId",how="left")
r["gap"]=r.Bandgap
print(f"band gap present for {r.gap.notna().sum()} of {len(r)} recomputed materials "
      f"(infinite values screened out)")

print("(a) flip rate by band-gap bin x 3d, and within bin by anion class")
bins=[-1e-9,1e-9,0.5,1.0,3.0,99]; labs=["metal (0)","0-0.5","0.5-1","1-3",">3"]
r["gb"]=pd.cut(r.gap,bins,labels=labs)
tab=r.pivot_table(index="gb",columns="d3",values="flip",aggfunc=["size","mean"],observed=True)
print(f"{'gap bin':12}{'n no-3d':>9}{'flip':>8}{'n 3d':>9}{'flip':>8}{'diff_pp':>9}")
rows=[]
for g in labs:
    sub=r[r.gb==g]
    a=sub[~sub.d3]; b=sub[sub.d3]
    if len(a)<100 or len(b)<100: continue
    rows.append(dict(gap_bin=g,n_no3d=len(a),flip_no3d=round(a.flip.mean()*100,1),
                     n_3d=len(b),flip_3d=round(b.flip.mean()*100,1),
                     diff_pp=round((b.flip.mean()-a.flip.mean())*100,1)))
    print(f"{g:12}{len(a):9d}{a.flip.mean()*100:7.1f}%{len(b):9d}{b.flip.mean()*100:7.1f}%"
          f"{(b.flip.mean()-a.flip.mean())*100:+9.1f}")
pd.DataFrame(rows).to_csv("results/bandgap_bins_3d.csv",index=False)
print("\n   wide-gap (>3 eV) materials by anion class, no-3d only:")
wg=r[(r.gb==">3")&(~r.d3)]
print(wg.groupby("anion").agg(n=("flip","size"),flip=("flip","mean")).query("n>=100")
      .assign(flip=lambda d:(d.flip*100).round(1)).to_string())

print("\n(b) directional shift on the PRIMARY cell definition (quaternary, no Ln)")
m=r.merge(s,on="MaterialId",how="inner"); m=m[m.pbe<=0].copy()
m["delta"]=(m.scan-m.pbe)*1000
q=m[(m.n_el==4)&(~m.Ln)]
print(f"{'cell':24}{'n':>7}{'mean':>8}{'median':>8}{'sd':>7}{'Q1':>8}{'Q3':>8}{'share>0':>9}")
rows=[]
for a in ("oxide","pnictide","intermetallic","boride/carbide/silicide","halide","chalcogenide"):
    for d in (True,False):
        sub=q[(q.anion==a)&(q.d3==d)]
        if len(sub)<200: continue
        lab=a[:16]+("+3d" if d else "")
        rows.append(dict(anion=a,d3=d,n=len(sub),mean=round(sub.delta.mean(),1),
            median=round(sub.delta.median(),1),sd=round(sub.delta.std(),1),
            q1=round(sub.delta.quantile(.25),1),q3=round(sub.delta.quantile(.75),1),
            share_up=round((sub.delta>0).mean()*100,1)))
        print(f"{lab:24}{len(sub):7d}{sub.delta.mean():+8.1f}{sub.delta.median():+8.1f}"
              f"{sub.delta.std():7.1f}{sub.delta.quantile(.25):+8.1f}{sub.delta.quantile(.75):+8.1f}"
              f"{(sub.delta>0).mean()*100:8.1f}%")
pd.DataFrame(rows).to_csv("results/shift_primary_cells.csv",index=False)

"""Referee round 4, item 4: the sharpest available test of the Hubbard-U reading.
MP2020 applies +U to V Cr Mn Fe Co Ni W Mo in OXIDES AND FLUORIDES only. Within the
halide class, fluorine-containing materials therefore receive +U on those metals and
Cl/Br/I-only materials do not. Same anion class, same elements, correction on or off.
If +U suppresses the disagreement where applied, the Fe excess must be SMALLER in the
F-containing subset.
Note that MP2020 anion corrections (O S F Cl Br I N H Se Si Sb Te) are per-anion
constants: they apply equally to Fe-containing and Fe-free materials of the same class
and therefore cancel in a within-anion excess. +U is the only correction that acts
differentially on the cation, which is what makes this contrast interpretable.
"""
import ast
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
UEL=set("V Cr Mn Fe Co Ni W Mo".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C])
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["flip"]=(r[C]>0).astype(int)
hal=r[r.anion=="halide"].copy()
hal["hasF"]=els[r.anion=="halide"].map(lambda e:"F" in e).values
print("Halides only. F present => +U applied to U-list metals; Cl/Br/I only => no +U.\n")
print(f"{'subset':16}{'metal':8}{'n_with':>8}{'flip_with':>11}{'flip_without':>14}{'excess_pp':>11}")
rows=[]
for lab,sub in (("F-containing",hal[hal.hasF]),("Cl/Br/I only",hal[~hal.hasF])):
    e_sub=els[sub.index]
    for sym in ("Fe","Co","Mn","any U-list"):
        if sym=="any U-list": has=e_sub.map(lambda e:bool(e&UEL)).values
        else: has=e_sub.map(lambda e,s=sym: s in e).values
        a,b=sub[has],sub[~has]
        if len(a)<60: 
            print(f"{lab:16}{sym:8}{len(a):8d}{'(n too small)':>36}"); continue
        exc=(a.flip.mean()-b.flip.mean())*100
        rows.append(dict(subset=lab,metal=sym,n_with=len(a),flip_with=round(a.flip.mean()*100,1),
                         flip_without=round(b.flip.mean()*100,1),excess_pp=round(exc,1)))
        print(f"{lab:16}{sym:8}{len(a):8d}{a.flip.mean()*100:10.1f}%{b.flip.mean()*100:13.1f}%{exc:+11.1f}")
pd.DataFrame(rows).to_csv("results/halide_u_split.csv",index=False)
print("\n-> results/halide_u_split.csv")
print("Prediction if +U suppresses: excess smaller in the F-containing subset.")

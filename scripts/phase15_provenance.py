"""Referee round 4, item: selection on PBE stability ("winner's curse").

GNoME candidates were filtered for PBE stability by a PBE-trained model, so the released
set is enriched for structures where PBE is anomalously favourable. That alone predicts
elevated flip rates relative to known materials, with no chemistry-functional interaction.
The propensity analysis in Section 6 models selection INTO the r2SCAN subset, not
selection into the stable set, so it does not address this.

Two probes the release supports:
 (1) generation pipeline, from the 'Data Directory' path. GNoME used symmetry-aware
     partial substitution and random structure search; the path segment marks the run,
     and 'wy' runs are the Wyckoff/symmetry-aware ones. Substituted structures sit near
     known chemistry; searched structures do not. Different exposure to the curse.
 (2) 'Is Train': whether the material was in the GNN's training data. Predictions on
     held-out chemistry carry more model-selection pressure than memorised ones.
If the chemistry gradient is a selection artefact it should differ sharply across these
splits. If it is functional disagreement it should not.
"""
import ast, re
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
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
s=pd.read_csv("data/stable_materials_summary.csv",
              usecols=["MaterialId","Data Directory","Is Train",C]).rename(columns={C:"pbe"})
s["pbe"]=pd.to_numeric(s.pbe,errors="coerce")
def run_of(p):
    m=re.search(r"/borg/([^/]+)/",str(p)); return m.group(1) if m else "unknown"
s["run"]=s["Data Directory"].map(run_of)
s["wy"]=s.run.str.contains("wy")
r=pd.read_csv("data/stable_materials_r2scan.csv",usecols=["MaterialId","Elements",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).rename(columns={C:"scan"})
els=r["Elements"].map(ast.literal_eval).map(set)
r["anion"]=els.map(ac); r["d3"]=els.map(lambda e:bool(e&TM3D)); r["n_el"]=els.map(len)
m=r.merge(s,on="MaterialId",how="inner"); m=m[m.pbe<=0].copy()
m["flip"]=(m.scan>0).astype(int); m["hull"]=-m.pbe
print(f"matched & PBE-stable: {len(m)}")
print("\npipeline split (Wyckoff/symmetry-aware vs other):")
print(m.groupby("wy").agg(n=("flip","size"),flip=("flip","mean"),
                          med_hull=("hull","median")).assign(flip=lambda d:(d.flip*100).round(2)).to_string())
print("\nIs Train split:")
print(m.groupby("Is Train").agg(n=("flip","size"),flip=("flip","mean"),
                          med_hull=("hull","median")).assign(flip=lambda d:(d.flip*100).round(2)).to_string())

print("\nHeadline contrast computed separately inside each split")
print("(if the gradient is a selection artefact it should not be stable across them):")
rows=[]
for lab,col in (("pipeline","wy"),("in training set","Is Train")):
    for val in (True,False):
        sub=m[m[col]==val]
        hi=sub[(sub.anion=="oxide")&sub.d3]; lo=sub[(sub.anion=="intermetallic")&~sub.d3]
        if len(hi)<150 or len(lo)<150: continue
        g=(hi.flip.mean()-lo.flip.mean())*100
        rows.append(dict(split=lab,value=str(val),n_hi=len(hi),hi=round(hi.flip.mean()*100,1),
                         n_lo=len(lo),lo=round(lo.flip.mean()*100,1),gap_pp=round(g,1)))
        print(f"  {lab:16}={str(val):5}  oxide+3d {hi.flip.mean()*100:5.1f}% (n={len(hi):6d})   "
              f"interm {lo.flip.mean()*100:4.1f}% (n={len(lo):6d})   gap {g:+5.1f}pp")
pd.DataFrame(rows).to_csv("results/provenance_split.csv",index=False)
print("\n-> results/provenance_split.csv")

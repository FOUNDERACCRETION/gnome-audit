"""Is the 3d/oxide effect really the functional, or just MP2020 +U corrections?
MP applies U only to Co Cr Fe Mn Mo Ni V W, and only in oxides/fluorides.
If the effect survives among 3d elements that never get U, it is chemistry.
If it vanishes, the headline is a corrections artefact."""
import ast
import numpy as np, pandas as pd
RNG = np.random.default_rng(20261009)
U_EL = set("Co Cr Fe Mn Mo Ni V W".split())
TM3D_NOU = set("Sc Ti Cu Zn".split())          # 3d metals MP never applies U to
COL = "Decomposition Energy Per Atom MP"

df = pd.read_csv("data/stable_materials_r2scan.csv", index_col=0)
els = df["Elements"].map(ast.literal_eval).map(set)
df["n_el"] = els.map(len)
df["has_O"] = els.map(lambda e: "O" in e)
df["has_F"] = els.map(lambda e: "F" in e)
df["U_el"] = els.map(lambda e: bool(e & U_EL))
df["noU_3d"] = els.map(lambda e: bool(e & TM3D_NOU) and not (e & U_EL))
df[COL] = pd.to_numeric(df[COL], errors="coerce")
df = df.dropna(subset=[COL]); df["flip"] = df[COL] > 0

def ci(v, nb=1000):
    n=len(v); p=float(np.mean(v)); d=RNG.binomial(n,p,nb)/n
    return p, float(np.percentile(d,2.5)), float(np.percentile(d,97.5))

print("=== Does U explain the 3d effect? (oxides/fluorides = where MP applies U) ===")
rows=[]
for n_el in (3,4):
    for ox in (True, False):
        sub = df[(df.n_el==n_el) & (df.has_O | df.has_F)==ox] if False else \
              df[(df.n_el==n_el) & ((df.has_O | df.has_F)==ox)]
        for label, m in [("U-element", sub.U_el),
                         ("non-U 3d (Sc/Ti/Cu/Zn)", sub.noU_3d),
                         ("no 3d at all", ~sub.U_el & ~sub.noU_3d)]:
            cell = sub[m]
            if len(cell) < 400: continue
            p,lo,hi = ci(cell.flip.values)
            rows.append(dict(n_el=n_el, oxide_or_fluoride=ox, group=label,
                             n=len(cell), flip=p, lo=lo, hi=hi))
r = pd.DataFrame(rows)
print(r.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
r.to_csv("results/u_confound_test.csv", index=False)

print("\n=== Same, outside oxides/fluorides (U is NOT applied there) ===")
print("If U-elements still flip more here, the effect is not the correction scheme.")

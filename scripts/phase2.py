"""Phase 2: chemistry-resolved PBE->r2SCAN flip rates for GNoME.
Decision rules pre-registered in DECISION_RULES.md (commit 82d49db) before this ran.
"""
import ast, json, sys
import numpy as np, pandas as pd

RNG = np.random.default_rng(20261009)
DATA = sys.argv[1] if len(sys.argv)>1 else "data/stable_materials_r2scan.csv"
OUT = sys.argv[2] if len(sys.argv)>2 else "results"

LN = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
AN = set("Ac Th Pa U Np Pu".split())
TM3D = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
HAL = set("F Cl Br I".split())
CHAL = set("S Se Te".split())
PNIC = set("N P As Sb Bi".split())
POT = set("Li Na Mg Ge Ga".split())           # newer PAW potentials vs MP standard
PLUSU = set("Co Cr Fe Mn Mo Ni V W".split())  # MP2020 +U elements

COL_MP = "Decomposition Energy Per Atom MP"
COL_OQMD = "Decomposition Energy Per Atom MP OQMD"

def anion_class(e):
    if "O" in e: return "oxide"
    if e & HAL: return "halide"
    if e & CHAL: return "chalcogenide"
    if e & PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e & {"B", "C", "Si"}: return "boride/carbide/silicide"
    return "intermetallic"

print("loading", DATA, flush=True)
df = pd.read_csv(DATA, index_col=0)
els = df["Elements"].map(ast.literal_eval).map(set)
df["n_elements"] = els.map(len)
df["anion_class"] = els.map(anion_class)
df["has_Ln"] = els.map(lambda e: bool(e & LN))
df["has_An"] = els.map(lambda e: bool(e & AN))
df["has_3d"] = els.map(lambda e: bool(e & TM3D))
df["has_newpot"] = els.map(lambda e: bool(e & POT))
df["has_U"] = els.map(lambda e: bool(e & PLUSU))
for c in (COL_MP, COL_OQMD):
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna(subset=[COL_MP]).copy()
df["flip_mp"] = df[COL_MP] > 0
df["flip_oqmd"] = df[COL_OQMD] > 0
print(f"rows={len(df)}  global flip MP={df.flip_mp.mean()*100:.2f}%  "
      f"OQMD={df.flip_oqmd.mean()*100:.2f}%", flush=True)

def ci(flips, n_boot=1000):
    """Bootstrap 95% CI for a proportion. Resampling Bernoulli data is exactly
    Binomial(n, p_hat)/n, so this is the nonparametric bootstrap, done fast."""
    n = len(flips)
    if n == 0: return (np.nan, np.nan, np.nan)
    p = float(np.mean(flips))
    draws = RNG.binomial(n, p, n_boot) / n
    return p, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))

# ---------- PRIMARY: disjoint cells, anion_class x has_Ln x has_3d, within n_elements
rows = []
for n_el, g in df.groupby("n_elements"):
    if len(g) < 500: continue
    for (ac, ln, d3), cell in g.groupby(["anion_class", "has_Ln", "has_3d"]):
        if len(cell) < 500: continue
        p, lo, hi = ci(cell.flip_mp.values)
        p2, lo2, hi2 = ci(cell.flip_oqmd.values)
        rows.append(dict(n_elements=n_el, anion_class=ac, has_Ln=ln, has_3d=d3,
                         n=len(cell), flip_mp=p, lo=lo, hi=hi,
                         flip_oqmd=p2, lo_oqmd=lo2, hi_oqmd=hi2))
strat = pd.DataFrame(rows).sort_values(["n_elements", "flip_mp"], ascending=[True, False])
strat.to_csv(f"{OUT}/stratified.csv", index=False)
print("\n=== PRIMARY: disjoint cells (n>=500) ===", flush=True)
print(strat.to_string(index=False, float_format=lambda x: f"{x:.4f}"), flush=True)

# spread per n_elements level
spreads = {}
for n_el, g in strat.groupby("n_elements"):
    if len(g) < 2: continue
    top, bot = g.iloc[0], g.iloc[-1]
    spreads[int(n_el)] = dict(
        spread_pp=float((top.flip_mp - bot.flip_mp) * 100),
        disjoint_ci=bool(top.lo > bot.hi),
        top=f"{top.anion_class}/Ln={top.has_Ln}/3d={top.has_3d}",
        top_rate=float(top.flip_mp), bot=f"{bot.anion_class}/Ln={bot.has_Ln}/3d={bot.has_3d}",
        bot_rate=float(bot.flip_mp), n_cells=int(len(g)))

# ---------- SECONDARY: Ln effect, controlling for anion_class and n_elements
ln_rows = []
for (n_el, ac, d3), g in df.groupby(["n_elements", "anion_class", "has_3d"]):
    a, b = g[g.has_Ln], g[~g.has_Ln]
    if len(a) < 500 or len(b) < 500: continue
    pa, loa, hia = ci(a.flip_mp.values); pb, lob, hib = ci(b.flip_mp.values)
    ln_rows.append(dict(n_elements=n_el, anion_class=ac, has_3d=d3,
                        n_Ln=len(a), rate_Ln=pa, lo_Ln=loa, hi_Ln=hia,
                        n_noLn=len(b), rate_noLn=pb, lo_noLn=lob, hi_noLn=hib,
                        delta_pp=(pa - pb) * 100, disjoint=bool(hia < lob or hib < loa)))
ln_df = pd.DataFrame(ln_rows)
if len(ln_df): ln_df.to_csv(f"{OUT}/lanthanide_effect.csv", index=False)
print("\n=== SECONDARY: lanthanide effect, matched cells ===", flush=True)
print(ln_df.to_string(index=False, float_format=lambda x: f"{x:.4f}") if len(ln_df) else "no matched cells", flush=True)

# ---------- ROBUSTNESS 1: threshold curves
tr = []
for t in np.arange(0, 0.051, 0.005):
    r = dict(threshold_eV=round(float(t), 3), overall=float((df[COL_MP] > t).mean()))
    for ac, g in df.groupby("anion_class"):
        if len(g) >= 500: r[ac] = float((g[COL_MP] > t).mean())
    tr.append(r)
thr = pd.DataFrame(tr)
thr.to_csv(f"{OUT}/threshold_curves.csv", index=False)
print("\n=== ROBUSTNESS: flip rate vs threshold ===", flush=True)
print(thr.to_string(index=False, float_format=lambda x: f"{x:.4f}"), flush=True)

# ---------- ROBUSTNESS 2: potential and +U confounds
conf = []
for name, mask in [("newpot_LiNaMgGeGa", df.has_newpot), ("plusU_elements", df.has_U)]:
    for val in (True, False):
        sub = df[mask == val]
        if len(sub) < 500: continue
        p, lo, hi = ci(sub.flip_mp.values)
        conf.append(dict(flag=name, value=val, n=len(sub), flip_mp=p, lo=lo, hi=hi))
conf_df = pd.DataFrame(conf); conf_df.to_csv(f"{OUT}/confounds.csv", index=False)
print("\n=== ROBUSTNESS: confound splits ===", flush=True)
print(conf_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"), flush=True)

# ---------- ROBUSTNESS 3: near-hull concentration (Note 2's claim)
bins = [-np.inf, -0.1, -0.05, -0.02, -0.01, 0.0]
df["pbe_bin"] = pd.cut(df[COL_MP], bins)  # NB: r2SCAN column; PBE bins added in Phase 3
nh = df.groupby("pbe_bin", observed=True).agg(n=("flip_mp", "size")).reset_index()
nh.to_csv(f"{OUT}/hull_distance_bins.csv", index=False)

# ---------- DECISION
best = max(spreads.items(), key=lambda kv: kv[1]["spread_pp"]) if spreads else None
verdict = "INCONCLUSIVE"
if best:
    s = best[1]
    if s["spread_pp"] >= 5 and s["disjoint_ci"]: verdict = "GO"
    elif s["spread_pp"] < 3: verdict = "FOLD"
    else: verdict = "GREY"
res = dict(rows=len(df), global_flip_mp=float(df.flip_mp.mean()),
           global_flip_oqmd=float(df.flip_oqmd.mean()), spreads=spreads, verdict=verdict)
with open(f"{OUT}/phase2_summary.json", "w") as f: json.dump(res, f, indent=2)
print("\n=== VERDICT:", verdict, "===", flush=True)
print(json.dumps(spreads, indent=2), flush=True)

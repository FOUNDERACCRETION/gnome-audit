#!/usr/bin/env python3
"""Regenerate every result in results/ and check the paper's headline numbers.

    python run_all.py                 run all phases, then check
    python run_all.py --check         only check existing results/ against the paper
    python run_all.py --skip a,b      skip named scripts (without .py)

Each script's stdout goes to results/<script>.log. Scripts that need the Materials
Project API are skipped unless MP_API_KEY (or MP_API_KEY_FILE) is set.
"""
import argparse, csv, json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)

ORDER = [
    "phase2", "phase2b", "phase3_selection", "phase5_budget", "phase6_extended",
    "phase6b_logit", "phase7_cluster", "phase8_referee", "phase9_selection",
    "phase10_robust", "phase11_hostile", "phase12_budget_baselines",
    "phase13_halide_split", "phase14_reference_offsets", "phase14b_additivity",
    "phase15_provenance", "phase16_excess_ci", "phase17_validation_deff",
    "phase18_gap_delta", "phase19_logit_ipw", "phase20_audit",
    "phase21_openshell", "phase22_populations",
    "phase24_budget_strong", "phase25_register", "phase29_reconcile",
    "phase31_stability", "phase31b_fe_stability", "phase32_bh_identity",
    "phase33_percentile_conv", "phase34_converged_intervals",
    "phase35_coarse_cluster",
    "q2_known_baseline", "q2b_completeness",
    "phase26_diag", "phase27_ame_ci", "phase28_standardise", "phase30_matched_basis",
    "phase37_matched_cluster", "phase38_resampling_equivalence",
    "phase40_allref_cells", "phase41_hull_continuous", "phase43_context_identify",
    "make_figures",
]
NEEDS_MP = {"q2_known_baseline", "q2b_completeness"}
# phase26_diag and phase28_standardise read the cached Materials Project pulls but make no
# API call; they run from the stored CSVs in CI.

# Every number in the abstract. (file, how to extract, paper value, tolerance)
def _json(path, *keys):
    v = json.load(open(path))
    for k in keys: v = v[k]
    return v

def _row(path, where, col):
    for r in csv.DictReader(open(path)):
        if all(str(r[k]) == str(v) for k, v in where.items()):
            return float(r[col])
    raise KeyError(f"{where} not in {path}")

CHECKS = [
    ("global flip, MP reference (%)",
     lambda: 100*_json("results/phase2_summary.json", "global_flip_mp"), 13.04, 0.1),
    ("global flip, MP+OQMD reference (%)",
     lambda: 100*_json("results/phase2_summary.json", "global_flip_oqmd"), 15.93, 0.1),
    ("quaternary spread, top-bottom cell (pp)",
     lambda: _json("results/phase2_summary.json", "spreads", "4", "spread_pp"), 26.7, 0.1),
    ("quaternary oxide/3d/no-Ln cell (%)",
     lambda: 100*_row("results/stratified.csv",
                      {"n_elements": 4, "anion_class": "oxide", "has_Ln": False, "has_3d": True}, "flip_mp"), 30.7, 0.1),
    ("quaternary intermetallic/no-3d/no-Ln cell (%)",
     lambda: 100*_row("results/stratified.csv",
                      {"n_elements": 4, "anion_class": "intermetallic", "has_Ln": False, "has_3d": False}, "flip_mp"), 4.0, 0.1),
    ("cluster-robust OR, has_U",
     lambda: _row("results/logit_cluster_robust.csv", {"feature": "has_U"}, "OR"), 3.45, 0.02),
    ("cluster-robust OR, has_U, lower",
     lambda: _row("results/logit_cluster_robust.csv", {"feature": "has_U"}, "lo_B"), 3.13, 0.02),
    ("cluster-robust OR, has_U, upper",
     lambda: _row("results/logit_cluster_robust.csv", {"feature": "has_U"}, "hi_B"), 3.80, 0.02),
    # Section 4.16: the ORDERING of the four policies is the claim, not any one value.
    ("budget: composition at 20% (%)",
     lambda: 100*_row("results/budget_baselines.csv", {"budget": 0.2}, "chemistry"), 37.6, 0.2),
    ("budget: hull distance at 20% (%)",
     lambda: 100*_row("results/budget_baselines.csv", {"budget": 0.2}, "hull"), 40.4, 0.2),
    ("budget: both at 20% (%)",
     lambda: 100*_row("results/budget_baselines.csv", {"budget": 0.2}, "combined"), 45.8, 0.2),
]

def _excess(el, anion):
    for r in csv.DictReader(open("results/within_anion_excess.csv")):
        if r["element"] == el and r["anion"] == anion:
            return float(r["excess_pp"])
    raise KeyError(f"{el}/{anion} not in within_anion_excess.csv")

# Section 4.3, tab:feexcess. The claim the whole section rests on is that the Fe excess is
# SMALLEST in oxides (where PBE gets a +U) and largest where it does not.
CHECKS += [
    ("Fe excess, boride/carbide/silicide (pp)",
     lambda: _excess("Fe", "boride/carbide/silicide"), 41.0, 0.1),
    ("Fe excess, pnictide (pp)", lambda: _excess("Fe", "pnictide"), 35.3, 0.1),
    ("Fe excess, intermetallic (pp)", lambda: _excess("Fe", "intermetallic"), 34.2, 0.1),
    ("Fe excess, oxide (pp, must be the smallest)", lambda: _excess("Fe", "oxide"), 5.0, 0.1),
    ("Co excess, boride/carbide/silicide (pp)",
     lambda: _excess("Co", "boride/carbide/silicide"), 34.0, 0.1),
    ("Co excess, oxide (pp)", lambda: _excess("Co", "oxide"), 4.6, 0.1),
    ("Ni excess, oxide (pp)", lambda: _excess("Ni", "oxide"), 10.3, 0.1),
    ("Ni excess, intermetallic (pp)", lambda: _excess("Ni", "intermetallic"), 2.0, 0.1),
]

def _ordering_holds():
    """The budget claim is an ordering: random < composition < hull < both, at every
    budget. A point check on one row would miss a reversal elsewhere."""
    for r in csv.DictReader(open("results/budget_baselines.csv")):
        v = [float(r["random"]), float(r["chemistry"]), float(r["hull"]), float(r["combined"])]
        if not (v[0] < v[1] < v[2] < v[3]):
            return 0.0
    return 1.0

def _excess_oxide_smallest():
    """Section 4.3 turns on the Fe excess being SMALLEST in oxides. Check the claim,
    not just the number."""
    rows = [r for r in csv.DictReader(open("results/within_anion_excess.csv"))
            if r["element"] == "Fe"]
    ox = [float(r["excess_pp"]) for r in rows if r["anion"] == "oxide"][0]
    return 1.0 if all(ox <= float(r["excess_pp"]) for r in rows) else 0.0

CHECKS += [
    ("[result] budget ordering random<composition<hull<both at every budget",
     _ordering_holds, 1.0, 0.0),
    ("[result] Fe excess is smallest in oxides of all anion classes",
     _excess_oxide_smallest, 1.0, 0.0),
    ("per-element contrasts surviving BH (of 83)",
     lambda: sum(r["sig_BH_05"] == "True"
                 for r in csv.DictReader(open("results/per_element_clustered.csv"))), 53, 0),
]

def known_aggregate():
    n = f = 0.0
    for r in csv.DictReader(open("results/known_baseline.csv")):
        n += float(r["n"]); f += float(r["n"]) * float(r["flip_known"])
    return 100 * f / n
def _allref(cell, field):
    import csv as _csv
    for r_ in _csv.DictReader(open("results/allref_cells.csv")):
        if r_["reference"] == "All" and r_["cell"] == cell:
            return float(r_[field])
    raise KeyError(cell)

CHECKS += [
    ("All-reference 4-el 3d oxide cell rate (%)",
     lambda: _allref("oxide, no Ln, 3d", "rate_pct"), 37.4, 0.1),
    ("All-reference 4-el intermetallic cell rate (%)",
     lambda: _allref("intermetallic, no Ln/3d", "rate_pct"), 14.5, 0.1),
    ("All-reference headline cell spread (pp)",
     lambda: _allref("spread (oxide - intermetallic, pp)", "rate_pct"), 22.9, 0.1),
    ("All-reference headline cell ratio",
     lambda: _allref("ratio (oxide / intermetallic)", "rate_pct"), 2.58, 0.02),
    ("[result] headline cell spread still excludes zero under the All reference",
     lambda: float(_allref("spread (oxide - intermetallic, pp)", "lo") > 0), 1.0, 0.0),
]

def _hullc(model, feature, field="odds_ratio"):
    import csv as _csv
    for r_ in _csv.DictReader(open("results/hull_continuous.csv")):
        if r_["model"] == model and r_["feature"] == feature:
            return float(r_[field])
    raise KeyError((model, feature))

_SPL = "natural spline in PBE decomposition energy"
CHECKS += [
    ("continuous-adjustment population (PBE-stable, matched)",
     lambda: _hullc("population", "n_matched_with_pbe_decomp"), 117018, 0.5),
    ("U-list OR, no hull-distance term",
     lambda: _hullc("no hull-distance term", "has_U"), 3.533, 0.01),
    ("U-list OR, natural spline in PBE depth",
     lambda: _hullc(_SPL, "has_U"), 3.431, 0.01),
    ("pseudo-R2 without hull distance",
     lambda: _hullc("no hull-distance term", "mcfadden_r2"), 0.0682, 0.001),
    ("pseudo-R2 with the hull-distance spline",
     lambda: _hullc(_SPL, "mcfadden_r2"), 0.1556, 0.001),
    ("[result] the U-list association survives continuous hull-distance adjustment",
     lambda: float(_hullc(_SPL, "has_U", "or_lo") > 1.0), 1.0, 0.0),
]

def _ctx(el, field):
    import csv as _csv
    for r_ in _csv.DictReader(open("results/context_spread_ci.csv")):
        if r_["element"] == el:
            return float(r_[field])
    raise KeyError(el)

def _jt(el, field):
    import csv as _csv
    for r_ in _csv.DictReader(open("results/context_joint_test.csv")):
        if r_["element"] == el:
            return float(r_[field])
    raise KeyError(el)

def _jtb(el):
    import csv as _csv
    for r_ in _csv.DictReader(open("results/context_joint_test.csv")):
        if r_["element"] == el:
            return r_["significant_bh_0.05"].strip().lower() == "true"
    raise KeyError(el)

def _ident(key):
    import csv as _csv
    for r_ in _csv.DictReader(open("results/context_identification.csv")):
        if r_["element"] == "(design)" and r_["context_a"] == key:
            return r_["context_b"]
    raise KeyError(key)

CHECKS += [
    ("[impl]   context model design is full rank (no minimum-norm choice in play)",
     lambda: float(_ident("rank") == "113 of 113"), 1.0, 0.0),
    ("[impl]   within-element contrasts invariant under reparameterisation (meV/atom)",
     lambda: float(_ident("max_reparam_contrast_diff_meV")), 0.0, 1e-6),
    ("Fe context spread (meV/atom)",   lambda: _ctx("Fe", "spread_mev"), 239, 0.6),
    ("Fe context spread, CI lower",    lambda: _ctx("Fe", "lo"), 161, 1.0),
    ("Fe context spread, CI upper",    lambda: _ctx("Fe", "hi"), 334, 1.0),
    ("Ni context spread (meV/atom)",   lambda: _ctx("Ni", "spread_mev"), 92, 0.6),
    ("Co context spread (meV/atom)",   lambda: _ctx("Co", "spread_mev"), 36, 0.6),
    ("context bootstrap replicates pooled over seeds",
     lambda: _ctx("Fe", "replicates"), 1965, 0.5),
    ("Ni joint-test chi2 (interactions equal)",  lambda: _jt("Ni", "wald_chi2"), 94.75, 0.02),
    ("Mn joint-test chi2 (interactions equal)",  lambda: _jt("Mn", "wald_chi2"), 71.63, 0.02),
    ("Fe joint-test chi2 (interactions equal)",  lambda: _jt("Fe", "wald_chi2"), 70.31, 0.02),
    ("Cr joint-test chi2 (interactions equal)",  lambda: _jt("Cr", "wald_chi2"), 31.71, 0.02),
    ("Co joint-test chi2 (interactions equal)",  lambda: _jt("Co", "wald_chi2"),  8.15, 0.02),
    ("[result] four of five elements reject equal interactions under BH",
     lambda: float(sum(_jtb(e) for e in ("Fe","Co","Ni","Mn","Cr"))), 4.0, 0.0),
    ("[result] cobalt does not reject equal interactions",
     lambda: float(not _jtb("Co")), 1.0, 0.0),
]

CHECKS_MP = [
    ("known-materials aggregate flip, cells with n>=100 (%)", known_aggregate, 14.88, 0.3),
    # --- added for the October 2026 revision register ---
    ("budget: eligible prototype cluster pool",
     lambda: _row("results/populations.csv", {"quantity": "budget_cluster_pool"}, "value"), 7822, 0.5),
    ("budget: estimation clusters per split",
     lambda: _row("results/populations.csv", {"quantity": "budget_clusters_per_split"}, "value"), 3911, 0.5),
    ("matched entries that were not PBE-stable",
     lambda: _row("results/populations.csv", {"quantity": "matched_not_pbe_stable"}, "value"), 25, 0.5),
    ("of those, inside a primary cell",
     lambda: _row("results/populations.csv", {"quantity": "of_those_in_primary_cells"}, "value"), 24, 0.5),
    ("flexible hull model minus closest-first ranking at 20% (pp)",
     lambda: _row("results/budget_strong_paired.csv", {"contrast": "hull_model-hull_rank"}, "mean_pp"), -0.06, 0.15),
    ("combined minus flexible hull model at 20% (pp)",
     lambda: _row("results/budget_strong_paired.csv", {"contrast": "combined-hull_model"}, "mean_pp"), 6.00, 0.3),
    ("cross-fitted combined minus flexible hull model at 20% (pp)",
     lambda: _row("results/budget_strong_paired.csv", {"contrast": "combined_cf-hull_model"}, "mean_pp"), 6.31, 0.3),
    ("[result] composition still adds after a nonparametric hull model (must be positive)",
     lambda: float(_row("results/budget_strong_paired.csv", {"contrast": "combined_cf-hull_model"}, "mean_pp") > 0), 1.0, 0.01),
    ("average marginal effect of U-list membership (pp)",
     lambda: _row("results/ame_ci.csv", {"term": "has_U"}, "ame_pp"), 15.4, 0.1),
    ("AME of U-list membership, lower bound (pp)",
     lambda: _row("results/ame_ci.csv", {"term": "has_U"}, "lo"), 14.3, 0.15),
    ("AME of U-list membership, upper bound (pp)",
     lambda: _row("results/ame_ci.csv", {"term": "has_U"}, "hi"), 16.6, 0.15),
    ("IPW weighted rate at the 1e-2 floor (%)",
     lambda: 100*_row("results/ipw_cells.csv", {"anion": "oxide", "d3b": "True"}, "ipw"), 27.2, 0.2),
    # --- added after the external review of 11 October 2026 ---
    ("smallest fitted propensity among weighted entries",
     lambda: _row("results/reconcile_checks.csv", {"quantity": "ps_min_eligible_recomputed"}, "value"), 1.519e-06, 1e-8),
    ("[impl]   untrimmed maximum weight = 1/that minimum",
     lambda: _row("results/reconcile_checks.csv", {"quantity": "untrimmed_max_weight"}, "value"), 658307.9, 1.0),
    ("eligible recomputed entries below the 1e-2 floor",
     lambda: _row("results/reconcile_checks.csv", {"quantity": "below_1e-2_floor"}, "value"), 26, 0.5),
    ("[impl]   hydride drop, Table 10 design: has_U odds ratio with hydrides",
     lambda: _row("results/reconcile_checks.csv", {"quantity": "t10_or_with_hydrides_has_U"}, "value"), 3.448, 0.002),
    ("[impl]   hydride drop, Table 10 design: has_3d_nonU with hydrides (must match Table 10)",
     lambda: _row("results/reconcile_checks.csv", {"quantity": "t10_or_with_hydrides_has_3d_nonU"}, "value"), 1.101, 0.002),
    ("hydride drop, Table 10 design: largest odds-ratio shift of any term",
     lambda: _row("results/reconcile_checks.csv", {"quantity": "t10_hydride_drop_max_shift"}, "value"), 0.023, 0.002),
    ("Fe oxide excess interval, lower endpoint at 100k replicates (pp)",
     lambda: _row("results/within_anion_excess_ci.csv", {"element": "Fe", "anion": "oxide"}, "lo"), 0.5, 0.051),
    ("Fe boride excess interval, lower endpoint (pp)",
     lambda: _row("results/within_anion_excess_ci.csv", {"element": "Fe", "anion": "boride/carbide/silicide"}, "lo"), 35.4, 0.051),
    ("per-element contrasts significant under ALL five seeds, at 20k replicates",
     lambda: _row("results/stability_checks.csv", {"quantity": "bh_core_all_seeds"}, "value"), 53, 0.5),
    ("[result] elements whose significance depends on the seed, at 20k replicates (must be none)",
     lambda: _row("results/stability_checks.csv", {"quantity": "bh_movers"}, "value"), 0, 0.5),
    ("percentile estimator converges: endpoint spread at 100k replicates (pp)",
     lambda: _row("results/percentile_convergence.csv", {"replicates": "100000"}, "lo_spread"), 0.035, 0.02),
    ("converged quaternary oxide+3d cell interval, lower (%)",
     lambda: 100*_row("results/cell_ci_converged.csv",
                      {"n_el": "4", "anion": "oxide", "Ln": "False", "d3": "True"}, "clust_lo"), 27.7, 0.051),
    ("converged quaternary oxide+3d cell interval, upper (%)",
     lambda: 100*_row("results/cell_ci_converged.csv",
                      {"n_el": "4", "anion": "oxide", "Ln": "False", "d3": "True"}, "clust_hi"), 33.8, 0.051),
    ("converged plain quaternary intermetallic interval, lower (%)",
     lambda: 100*_row("results/cell_ci_converged.csv",
                      {"n_el": "4", "anion": "intermetallic", "Ln": "False", "d3": "False"}, "clust_lo"), 3.3, 0.051),
    ("[result] headline cell intervals remain disjoint at the converged count",
     lambda: float(_row("results/cell_ci_converged.csv",
                        {"n_el": "4", "anion": "oxide", "Ln": "False", "d3": "True"}, "clust_lo")
                   > _row("results/cell_ci_converged.csv",
                        {"n_el": "4", "anion": "intermetallic", "Ln": "False", "d3": "False"}, "clust_hi")), 1.0, 0.01),
    ("converged Fe per-element interval, lower (%)",
     lambda: 100*_row("results/per_element_clustered.csv", {"element": "Fe"}, "clust_lo"), 39.8, 0.051),
    ("converged Fe per-element interval, upper (%)",
     lambda: 100*_row("results/per_element_clustered.csv", {"element": "Fe"}, "clust_hi"), 46.0, 0.051),
    ("matched cells with intervals excluding zero (cluster bootstrap)",
     lambda: float(sum(1 for r_ in csv.DictReader(open("results/known_cells_clustered_release.csv"))
                       if r_["excl0"] == "True")), 10, 0.5),
    ("matched cells with intervals excluding zero (material-level, superseded)",
     lambda: float(sum(1 for r_ in csv.DictReader(open("results/known_cells_clustered_release.csv"))
                       if r_["excl0_binom"] == "True")), 11, 0.5),
    ("matched cells where the cluster interval is wider than the independence-based one",
     lambda: float(sum(1 for r_ in csv.DictReader(open("results/known_cells_clustered_release.csv"))
                       if float(r_["widening"]) > 1.0)), 20, 0.5),
    ("[impl]   optimised cluster resampling agrees with an explicit row-gather resample",
     lambda: _row("results/resampling_equivalence.csv",
                  {"quantity": "equivalence_checks_failed"}, "value"), 0, 0.5),
    ("4-el 3d oxide matched difference, cluster CI lower (pp)",
     lambda: _row("results/known_cells_clustered_release.csv",
                  {"n_elements": "4", "anion_class": "oxide", "has_Ln": "False", "has_3d": "True"}, "lo"), 11.4, 0.051),
    ("[result] oxide and boride Fe intervals do not overlap",
     lambda: float(_row("results/within_anion_excess_ci.csv", {"element": "Fe", "anion": "boride/carbide/silicide"}, "lo")
                   > _row("results/within_anion_excess_ci.csv", {"element": "Fe", "anion": "oxide"}, "hi")), 1.0, 0.01),
]
# Set when the Materials Project phases did not run in this invocation, so that any check
# reading their output is reported as read from the stored file rather than regenerated.
MP_STORED = False

RUNTIMES = []

# phase43's spread intervals are a resumable accumulation: the paper's 1,965 replicates are
# 965 draws at seed 20261009 followed by 1,000 at seed 20261012 (the original first run was
# cut by its wall-clock budget at 965). Both invocations are fixed here, with the budget
# disabled, so the pooled draws are reproduced exactly rather than depending on run time.
MULTI = {"phase43_context_identify": [dict(RESET="1", SEED="20261009", NB="965", BUDGET="1e9"),
                                      dict(SEED="20261012", NB="1000", BUDGET="1e9")]}

def run(name):
    log = ROOT / "results" / f"{name}.log"
    t0 = time.time()
    with open(log, "w") as fh:
        for env_extra in MULTI.get(name, [{}]):
            env = dict(os.environ, **env_extra)
            p = subprocess.run([sys.executable, f"scripts/{name}.py"], stdout=fh, stderr=subprocess.STDOUT, env=env)
            if p.returncode != 0:
                break
    dt = time.time() - t0
    status = "ok" if p.returncode == 0 else f"FAILED (exit {p.returncode}, see {log})"
    print(f"  {name:22s} {dt:7.1f}s  {status}", flush=True)
    RUNTIMES.append(dict(phase=name, seconds=round(dt, 1), ok=p.returncode == 0))
    return p.returncode == 0

def write_runtimes():
    """Runtime is a reproducibility fact, so record it as a regenerated artefact rather
    than as a prose claim in the README that goes stale when a phase is added."""
    if not RUNTIMES:
        return
    import platform
    out = ROOT / "results" / "runtimes.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["phase", "seconds", "ok"])
        for r_ in RUNTIMES:
            w.writerow([r_["phase"], r_["seconds"], r_["ok"]])
        w.writerow(["TOTAL", round(sum(r_["seconds"] for r_ in RUNTIMES), 1), ""])
        w.writerow(["# platform", platform.platform(), ""])
        w.writerow(["# python", platform.python_version(), ""])
    print(f"\n  total {sum(r_['seconds'] for r_ in RUNTIMES)/60:.1f} min over "
          f"{len(RUNTIMES)} phases  ->  results/runtimes.csv", flush=True)

def check(mp=True):
    print("\nHeadline numbers vs paper.  [impl] = a property of the implementation that any"
          " legitimate run must satisfy;\n  [result] = an asserted finding, to be re-examined"
          " rather than forced if it fails;\n  everything else = a comparison against the"
          " recorded reference run, which is a staleness check and is\n  expected to be"
          " updated when a number legitimately changes.")
    bad = 0
    for label, fn, paper, tol in CHECKS + ([(l + "  [stored, not regenerated here]"
                                             if MP_STORED else l, f, p_, t_)
                                            for l, f, p_, t_ in CHECKS_MP] if mp else []):
        try:
            v = fn(); ok = abs(v - paper) <= tol
        except Exception as e:
            v, ok = float("nan"), False
            print(f"  ?  {label}: {e}")
            bad += 1; continue
        bad += (not ok)
        print(f"  {'ok' if ok else 'XX'} {label}: got {v:.3f}, paper {paper}")
    print("\nALL HEADLINE NUMBERS REPRODUCED" if not bad else f"\n{bad} MISMATCH(ES)")
    return bad == 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--skip", default="")
    a = ap.parse_args()
    skip = set(filter(None, a.skip.split(",")))
    have_key = bool(os.environ.get("MP_API_KEY") or os.environ.get("MP_API_KEY_FILE"))
    if not a.check:
        for d in ("data/stable_materials_r2scan.csv", "data/stable_materials_summary.csv"):
            if not Path(d).exists():
                sys.exit(f"missing {d}; see REPRODUCE.md section 2")
        Path("results").mkdir(exist_ok=True)
        print("Running phases:")
        failed = []
        for name in ORDER:
            if name in skip or (name in NEEDS_MP and not have_key):
                why = "--skip" if name in skip else "no MP_API_KEY"
                print(f"  {name:22s}          skipped ({why})")
                # Either way, if a Materials Project phase did not run, anything that
                # reads its output is a stored value and must be labelled as one.
                if name in NEEDS_MP:
                    MP_STORED = True
                continue
            if not run(name):
                failed.append(name)
        write_runtimes()
        if failed:
            # Do NOT let the check paper over a crashed phase: results/ may be committed,
            # so the check can pass on stale files while the code that made them is broken.
            print(f"\n{len(failed)} PHASE(S) FAILED: {', '.join(failed)}")
            check(mp=have_key or Path("results/known_baseline.csv").exists())
            sys.exit(1)
    ok = check(mp=have_key or Path("results/known_baseline.csv").exists())
    # the paper's tables must also match what the scripts just produced
    tex = Path("paper/gnome_audit_paper.tex")
    if tex.exists():
        print()
        r = subprocess.run([sys.executable, "scripts/check_paper.py", str(tex)])
        ok = ok and r.returncode == 0
        # check_paper.py only sees the tables. This catches a superseded number left in
        # prose, which is the failure that has recurred most often in this paper.
        print()
        r = subprocess.run([sys.executable, "scripts/check_stale.py", str(tex)])
        ok = ok and r.returncode == 0
    sys.exit(0 if ok else 1)

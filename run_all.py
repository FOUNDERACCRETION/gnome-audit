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
    "phase10_robust", "phase11_hostile",
    "q2_known_baseline", "q2b_completeness",
    "figures", "figures2", "figures3",
]
NEEDS_MP = {"q2_known_baseline", "q2b_completeness"}

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
    ("budget: chemistry-ranked catch at 20% (%)",
     lambda: 100*_row("results/budget_repeated.csv", {"budget": 0.2}, "chem_mean"), 37.1, 0.1),
    ("budget: chemistry-ranked sd at 20% (pp)",
     lambda: 100*_row("results/budget_repeated.csv", {"budget": 0.2}, "chem_sd"), 0.9, 0.1),
    ("budget: random catch at 20% (%)",
     lambda: 100*_row("results/budget_repeated.csv", {"budget": 0.2}, "random_mean"), 20.1, 0.1),
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

def known_aggregate():
    n = f = 0.0
    for r in csv.DictReader(open("results/known_baseline.csv")):
        n += float(r["n"]); f += float(r["n"]) * float(r["flip_known"])
    return 100 * f / n
CHECKS_MP = [
    ("known-materials aggregate flip, cells with n>=100 (%)", known_aggregate, 14.88, 0.3),
]
# Set when the Materials Project phases did not run in this invocation, so that any check
# reading their output is reported as read from the stored file rather than regenerated.
MP_STORED = False

def run(name):
    log = ROOT / "results" / f"{name}.log"
    t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.run([sys.executable, f"scripts/{name}.py"], stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    status = "ok" if p.returncode == 0 else f"FAILED (exit {p.returncode}, see {log})"
    print(f"  {name:22s} {dt:7.1f}s  {status}", flush=True)
    return p.returncode == 0

def check(mp=True):
    print("\nHeadline numbers vs paper:")
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
            if name in skip:
                print(f"  {name:22s}          skipped (--skip)"); continue
            if name in NEEDS_MP and not have_key:
                print(f"  {name:22s}          skipped (no MP_API_KEY)")
                MP_STORED = True
                continue
            if not run(name):
                failed.append(name)
        if failed:
            # Do NOT let the check paper over a crashed phase: results/ may be committed,
            # so the check can pass on stale files while the code that made them is broken.
            print(f"\n{len(failed)} PHASE(S) FAILED: {', '.join(failed)}")
            check(mp=have_key or Path("results/known_baseline.csv").exists())
            sys.exit(1)
    ok = check(mp=have_key or Path("results/known_baseline.csv").exists())
    sys.exit(0 if ok else 1)

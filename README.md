# gnome-audit

Code and derived results for:

> **Chemistry-resolved disagreement between PBE and r²SCAN stability verdicts in an
> AI-generated materials database.** Karthikeya (Accretion Labs), 2026. Preprint.

The paper asks one question of the GNoME r²SCAN recomputation set (Merchant et al.,
*Nature* 2023): does the rate at which a PBE "stable" verdict fails under r²SCAN depend
on chemistry? It does, by a factor of roughly seven between disjoint chemistry cells at
fixed element count, and the dependence survives clustering by structural prototype,
reweighting for recomputation selection, stratification by PBE hull distance, and a
comparison against known Materials Project materials. The paper then turns the
dependence into a verification-budget rule.

Everything in the paper is reproducible from two public CSV files and (for one section)
a free Materials Project API key. See [REPRODUCE.md](REPRODUCE.md).

## Layout

```
DECISION_RULES.md     pre-registered decision rules, committed before analysis (commit 82d49db)
DECISION.md           dated decision records written as each phase completed
scripts/              one script per analysis phase; see table below
results/              every CSV/PNG the paper draws a number or figure from
run_all.py            regenerates results/ end to end and checks headline numbers
paper/                LaTeX source and PDF of the current draft
.github/workflows/    continuous reproduction on a fresh runner (see below)
review/               adversarial referee report generated during drafting (round 3); the paper's
                      Sections 4.3-4.5 and 4.14 are the response to it
requirements.txt      pinned versions
data/                 NOT in the repository; two GNoME CSVs, see REPRODUCE.md
```

## Scripts, in the order `run_all.py` executes them

| Script | Paper section | What it does | Writes |
|---|---|---|---|
| `phase2.py` | §4.1, App. A, C | Disjoint cells (n_elements × anion × Ln × 3d), naive bootstrap CIs, threshold curves, GO/FOLD rule | `stratified.csv`, `lanthanide_effect.csv`, `threshold_curves.csv`, `confounds.csv`, `hull_distance_bins.csv`, `phase2_summary.json` |
| `phase2b.py` | §4.3 | Hubbard-*U* natural experiment (U-list vs non-U 3d vs no-3d, inside and outside O/F) | `u_confound_test.csv` |
| `phase3_selection.py` | §4.7, §5 | Recomputation coverage by cell and hull distance; post-stratification reweighting | `selection_coverage.csv`, `selection_by_hull_distance.csv`, `reweighted_cells.csv` |
| `phase5_budget.py` | §4.12 | Budget model, single fit/evaluate split | `budget_model.csv` |
| `phase6_extended.py` | §4.2, 4.6, 4.8, App. B | Per-element rates, band-gap bins, crystal system, cell size | `per_element.csv`, `bandgap_*.csv`, `crystal_system.csv`, `size_effect.csv`, `decomp_distribution.csv` |
| `phase6b_logit.py` | §4.9 | Logistic regression, IRLS, Fisher-information SEs | `logistic_coefficients.csv` |
| `phase7_cluster.py` | §3.3, 4.1 | Prototype clustering (space group × anonymised stoichiometry), ICC, cluster bootstrap | `cluster_bootstrap.csv` |
| `phase8_referee.py` | §4.2, 4.12 | Benjamini–Hochberg over elements; budget over 20 cluster-wise splits; per-element cluster CIs | `per_element_clustered.csv`, `budget_repeated.csv` |
| `phase9_selection.py` | §4.7 | Propensity model for recomputation; inverse-probability weighting | `ipw_cells.csv` |
| `phase10_robust.py` | §3.5, 4.9 | Cluster-robust sandwich SEs (three clusterings); alternative cluster definition | `logit_cluster_robust.csv` |
| `phase11_hostile.py` | §4.3, 4.4, 4.5, 4.14 | Pre-registered gate 6 (hull-distance strata); directional shift Δ; "All" reference; Fe/Co excess by anion; PBE-stable check | printed tables → `phase11.log` |
| `q2_known_baseline.py` | §4.11 | **Needs MP API key.** Pulls MP materials with GGA/GGA+U and r²SCAN thermo; same cells | `known_baseline.csv` |
| `q2b_completeness.py` | §4.11, 4.14, §7 | **Needs MP API key.** MP-side r²SCAN hull completeness; gradient within completeness terciles; hull-distance-matched Q2 | `mp_completeness.csv`, `mp_*_raw.csv` |
| `figures.py`, `figures2.py`, `figures3.py` | Figs 1–6 | Plots from `results/` only | `fig*.png` |

All random draws use `numpy.random.default_rng(20261009)`.

## The two things most likely to trip a re-run

1. The Materials Project `thermo_types` value is `r2SCAN`, case-sensitive. `R2SCAN`
   returns zero rows with no error.
2. `stable_materials_r2scan.csv` has `MaterialId` as its first column; reading it with
   `index_col=0` and `usecols` together silently drops it.

## Continuous reproduction

`.github/workflows/reproduce.yml` runs on every push to `master` and on demand: a fresh
Linux runner clones the repository, downloads the two GNoME CSVs from DeepMind's bucket,
regenerates `results/` from scratch (all phases except the two Materials Project pulls)
and runs `run_all.py --check`. The Actions tab holds the log; a green run means every
number in the abstract was regenerated on hardware the author does not own.

## Citing

Repository: https://github.com/FOUNDERACCRETION/gnome-audit. A `CITATION.cff` will be added when the preprint has an identifier.

## Licence

Code: MIT (see `LICENSE`). The GNoME data are released by Google DeepMind under
CC BY 4.0; Materials Project data under CC BY 4.0. Neither is redistributed here.

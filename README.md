# gnome-audit

Code and derived results for:

> **Chemistry-resolved disagreement between Materials Project-compatible PBE(+U) and
> r²SCAN stability verdicts in the GNoME materials database.**
> Karthikeya Aduru (Accretion Labs), 2026. Preprint; arXiv identifier to follow.

The paper asks one question of the GNoME r²SCAN recomputation set (Merchant et al.,
*Nature* 2023): does the rate at which a released PBE-level "stable" label fails under
r²SCAN depend on chemistry? It does — 4.0% to 30.7% between disjoint composition cells at
fixed element count, a factor of about seven — and the dependence survives clustering by
structural prototype, reweighting for recomputation selection, stratification and
continuous adjustment for PBE hull distance, two reference sets, and a matched comparison
against previously catalogued Materials Project materials. The paper then turns the
dependence into a verification-budget rule. It measures disagreement between two
workflows, not accuracy.

Everything in the paper is reproducible from two public CSV files and, for two phases, a
free Materials Project API key. **[REPRODUCE.md](REPRODUCE.md)** says how to regenerate
every number; **[VALIDATION.md](VALIDATION.md)** says what each automated check asserts
and what none of them can establish; **[REVIEWER.md](REVIEWER.md)** is for someone
checking whether the analysis is *correct*, which is a different job from checking that it
reproduces. Every defect found in this work so far was found by reading code against the
manuscript, not by running it, and the automated checks were green through all of them.

## Layout

```
DECISION_RULES.md     pre-committed decision rules, written before the stratified analysis (commit 82d49db)
DECISION.md           dated decision records written as each phase completed
scripts/              one script per analysis phase (table below), plus the three checkers
results/              every CSV, TeX fragment and figure the paper draws a number from
run_all.py            regenerates results/ end to end, then checks 84 headline quantities,
                      655 printed values (scripts/check_paper.py) and 25 superseded
                      literals (scripts/check_stale.py) against paper/gnome_audit_paper.tex
paper/                LaTeX source and PDF of the current draft; a split main+supplement
                      variant (_d32) is kept alongside
review/               adversarial referee report generated during drafting and the
                      author's responses; Sections 4.3–4.5 and 4.14 answer it
VALIDATION.md         what the checks assert; determinism notes; the live MP re-pull
REVIEWER.md           independent-review packet: where to look, what to doubt
REPRODUCE.md          step-by-step reproduction on Linux, macOS and Windows
CITATION.cff          how to cite the software
requirements.txt      pinned versions (pandas 2.3.3, numpy 2.2.6, matplotlib 3.10.9)
data/                 NOT in the repository; the two GNoME CSVs, see REPRODUCE.md
```

## Scripts, in the order `run_all.py` executes them

| Script | Paper | What it does | Writes |
|---|---|---|---|
| `phase2.py` | §4.1, App. A | Disjoint cells (element count × anion class × Ln × 3d), naive bootstrap intervals, threshold curves, GO/FOLD rule | `stratified.csv`, `lanthanide_effect.csv`, `threshold_curves.csv`, `confounds.csv`, `hull_distance_bins.csv`, `phase2_summary.json` |
| `phase2b.py` | §4.3 | Hubbard-*U* comparisons: U-list vs other-3d vs no-3d, inside and outside O/F chemistry | `u_confound_test.csv` |
| `phase3_selection.py` | §4.8, §5 | Recomputation coverage by cell and hull distance; post-stratification reweighting | `selection_coverage.csv`, `selection_by_hull_distance.csv`, `reweighted_cells.csv` |
| `phase5_budget.py` | §4.12 | Budget model, single fit/evaluate split | `budget_model.csv` |
| `phase6_extended.py` | §4.2, 4.7, 4.8 | Per-element rates, band-gap bins, crystal system, cell size | `per_element.csv`, `bandgap_*.csv`, `crystal_system.csv`, `size_effect.csv`, `decomp_distribution.csv` |
| `phase6b_logit.py` | §4.9 | Logistic regression, IRLS, Fisher-information SEs | `logistic_coefficients.csv` |
| `phase7_cluster.py` | §3.5 | Prototype clustering (space group × anonymised stoichiometry), ICC, cluster bootstrap | `cluster_bootstrap.csv` |
| `phase8_referee.py` | §4.2, 4.12 | Benjamini–Hochberg over elements; budget over 20 cluster-wise splits | `per_element_clustered.csv` (superseded by phase34), `budget_repeated.csv` |
| `phase9_selection.py` | §4.8 | Propensity model for recomputation; inverse-probability weighting | `ipw_cells.csv` |
| `phase10_robust.py` | §3.5, 4.9 | Cluster-robust sandwich SEs (three clusterings); alternative cluster definition | `logit_cluster_robust.csv` |
| `phase11_hostile.py` | §4.3, 4.5, 4.6, 4.14 | Hull-distance strata; directional shift Δ; "All" reference; within-anion Fe/Co/Ni excess | `within_anion_excess.csv`, tables → `phase11_hostile.log` |
| `phase12_budget_baselines.py` | §4.12 | Budget policies: random, composition, hull distance, both; 20 cluster-wise splits | `budget_baselines.csv` |
| `phase13_halide_split.py` | §4.3 | Fluorides (U applied) vs Cl/Br/I (not) within the halide class | `halide_u_split.csv` |
| `phase14_reference_offsets.py` | §4.4 | Additive per-element model of Δ and its collinearity with the cells | `reference_offsets.csv` |
| `phase14b_additivity.py` | §4.4 | Element × anion-class interactions in the fitted model | `context_dependence.csv` |
| `phase15_provenance.py` | §4.13 | Gradient across generation pipeline and training-set membership | `provenance_split.csv` |
| `phase16_excess_ci.py` | §4.3 | Cluster bootstrap intervals on the within-anion excess table (100,000 replicates) | `within_anion_excess_ci.csv` |
| `phase17_validation_deff.py` | §3.5, 4.1 | Element-count-matched validation against the source paper; design effect; drop-largest-clusters check | `phase17_validation_deff.log` |
| `phase18_gap_delta.py` | §4.6, 4.7 | Band gap in bins × 3d × anion; directional shift on the primary cells | `bandgap_bins_3d.csv`, `shift_primary_cells.csv` |
| `phase19_logit_ipw.py` | §4.9, 4.10 | Logistic model with anion interactions; likelihood-ratio test | `logit_interactions.csv` |
| `phase20_audit.py` | §3.7, 4.3, 4.7 | Verification of challenged numbers: band-gap accounting, lanthanide comparator, fluoride baselines | `phase20_audit.log` |
| `phase21_openshell.py` | §4.3 | U-list vs 3d-subset vs open-shell indicators, compared by likelihood | `phase21_openshell.log` |
| `phase22_populations.py` | §3.3, §5 | Population accounting; censoring bounds on the eligible release population | `phase22_populations.log` |
| `phase24_budget_strong.py` | §4.12, App. C | Budget test against a *modelled* hull distance; cross-fitted combined model | `budget_strong.csv`, `budget_strong_paired.csv` |
| `phase25_register.py` | §3.3, App. C | Population table and threshold curves at the register's resolution | `populations.csv`, `threshold_curves.csv` |
| `phase29_reconcile.py` | §4.8, App. C | Exact propensity minima, untrimmed weights, hydride-drop sensitivity | `reconcile_checks.csv` |
| `phase31_stability.py` | App. C | Seed stability of headline interval endpoints at 600 replicates | `phase31_stability.log` |
| `phase31b_fe_stability.py` | App. C | Seed stability of the Fe-in-oxides excess interval | `fe_oxide_seed_stability.csv` |
| `phase32_bh_identity.py` | §4.2, App. C | Which elements cross the BH threshold under five seeds at 20,000 replicates | `stability_checks.csv` |
| `phase33_percentile_conv.py` | App. C | Convergence of the percentile estimator with replicate count (400 → 100,000) | `percentile_convergence.csv` |
| `phase34_converged_intervals.py` | §4.1, 4.2, App. A, B | Every cell and per-element interval at the converged count (20,000) | `cell_ci_converged.csv`, `per_element_clustered.csv`, `cluster_bootstrap.csv`, `appendixA.tex` |
| `phase35_coarse_cluster.py` | §3.5 | Coarser clustering (stoichiometry only) at the converged count | `phase35_coarse_cluster.log` |
| `q2_known_baseline.py` | §4.11 | **Needs MP API key.** Pulls Materials Project materials with GGA/GGA+U and r²SCAN thermo records | `known_baseline.csv` |
| `q2b_completeness.py` | §4.11, 4.14 | **Needs MP API key.** MP-side r²SCAN hull completeness by anion class | `mp_completeness.csv`, `mp_gga_raw.csv`, `mp_r2scan_raw.csv` |
| `phase26_diag.py` | §4.8, 4.11 | IPW diagnostics; direct standardisation of the GNoME/MP comparison | `known_all_cells.csv` |
| `phase27_ame_ci.py` | §4.9 | Average marginal effects of the interaction model, cluster-bootstrap intervals (1,000 refits) | `ame_ci.csv` |
| `phase28_standardise.py` | §4.11, App. D | Standardised GNoME-minus-known differences over the matched cells | `known_all_cells.csv` |
| `phase30_matched_basis.py` | §4.11 | Matched comparison on the confirmed-transition basis | `known_baseline.csv` (columns) |
| `phase37_matched_cluster.py` | §4.11, App. D | Cluster bootstrap for the matched cells (GNoME side resampled by prototype) | `known_cells_clustered_release.csv`, `known_tab11.tex`, `known_all_cells.tex` |
| `phase38_resampling_equivalence.py` | §3.5 | Optimised cluster resampling agrees with an explicit row-gather resample | `resampling_equivalence.csv` |
| `phase40_allref_cells.py` | §4.14 | Headline cells under the "All" reference set, shared cluster draw | `allref_cells.csv` |
| `phase41_hull_continuous.py` | §4.5 | Composition model with a natural cubic spline in continuous PBE hull distance | `hull_continuous.csv` |
| `phase43_context_identify.py` | §4.4, App. C | Context model: identification, joint tests, cluster bootstrap over two seeds (965 + 1,000 draws) | `context_identification.csv`, `context_joint_test.csv`, `context_spread_ci.csv`, `context_spread_draws.npy` |
| `make_figures.py` | Figs 1–3 | The four figure PNGs, drawn from `results/` only, at print size | `fig1_…png` … `fig4_…png` |
| `check_paper.py` | all tables | Parses every number out of the LaTeX tables and prose and diffs it against `results/`: 655 values | exit status |
| `check_stale.py` | all text | Greps the manuscript for 25 superseded literals from earlier drafts | exit status |

`run_all.py` runs `phase43_context_identify.py` twice with fixed seeds (20261009 for 965
draws, then 20261012 for 1,000) to reproduce the paper's pooled 1,965 replicates. Seeds are
listed in the paper's Appendix C. `figures*.py`, `phase23_register.py`,
`phase36_ame_stability.py` and `phase42_pull_hull_depth.py` are retained for the record
but are not part of the driver and nothing in the paper depends on them.

## The three things most likely to trip a re-run

1. The Materials Project `thermo_types` value is `r2SCAN`, case-sensitive. `R2SCAN`
   returns zero rows with no error.
2. `stable_materials_r2scan.csv` has `MaterialId` as its first column; reading it with
   `index_col=0` and `usecols` together silently drops it.
3. numpy 2.2.6 has no wheel for Python 3.14; use Python 3.12 or 3.13. On Windows,
   `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` before activating the venv.

## Continuous reproduction

`.github/workflows/reproduce.yml` runs on every push to `master` and on demand: a fresh
Linux runner clones the repository, downloads the two GNoME CSVs from DeepMind's bucket,
deletes the committed contents of `results/` (keeping only the two Materials Project
pulls), regenerates everything, and runs the three checks. A green run means every
headline number and every value printed in the paper's tables and prose was regenerated on
hardware the author does not own and matched. The second check is the one that catches
transcription drift between an analysis and the manuscript describing it; it has found a
real off-by-one in the hull-completeness table and a replicate-count mismatch in the
context-model intervals. The two Materials Project analyses need an API key and do not run
there; their stored output is kept, and any check that reads it is labelled
`[stored, not regenerated here]` in the log. With a key set (`MP_API_KEY_FILE`), the full
pipeline has reproduced every headline number on four machines; the Materials Project
changes over time, and VALIDATION.md records the drift observed on a live re-pull.

## Citing

Until the preprint has an identifier, cite the repository:
https://github.com/FOUNDERACCRETION/gnome-audit (see `CITATION.cff`).

## Licence

Code: MIT (see `LICENSE`). The GNoME data are released by Google DeepMind under
CC BY-NC 4.0 and are not redistributed here. The cached Materials Project pulls in
`results/` are derived from Materials Project data (CC BY 4.0); see
`results/PROVENANCE.md`.

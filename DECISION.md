# Phase 2 decision record — GNoME PBE/r2SCAN chemistry-resolved flip rates
Run 2026-10-09, against rules pre-registered in DECISION_RULES.md (commit 82d49db).
Data: stable_materials_r2scan.csv, 139,541 rows. Global flip rate vs MP hull 13.04%,
vs MP+OQMD 15.93%. Reproduces published 86.8% quaternary survival.

## VERDICT: GO (standalone paper)

Rule required >=5pp spread with non-overlapping bootstrap CIs within a fixed
n_elements level. Observed:

- n_elements=4: spread 26.7pp across 23 disjoint cells, CIs disjoint.
  highest  oxide + 3d, no Ln        30.7%  [28.97, 32.57]  n=2,478
  lowest   intermetallic, no Ln/3d   4.0%  [ 3.48,  4.54]  n=5,025
- n_elements=3: spread 23.9pp across 17 cells, CIs disjoint.

Spread is 5x the GO threshold. Decision is not marginal.

## Secondary claim (lanthanides): PARTIALLY REFUTED — and more interesting for it

Pre-registered claim was that Ln materials flip at or below the global rate.
True in most matched cells, but it reverses by anion class:
  halide  n=4, with 3d   Ln 9.1% vs non-Ln 22.4%   -13.3pp  (Ln protective)
  oxide   n=4, with 3d   Ln 18.1% vs non-Ln 30.7%  -12.6pp  (Ln protective)
  pnictide n=4, with 3d  Ln 30.1% vs non-Ln 17.8%  +12.3pp  (Ln harmful)
  pnictide n=4, no 3d    Ln 13.9% vs non-Ln  6.0%   +7.9pp  (Ln harmful)
Sign flips with anion. Report as stated; do not collapse to one Ln number.

## Critical confound test: the effect is NOT the MP2020 +U correction

MP applies U only to Co Cr Fe Mn Mo Ni V W, and only in oxides/fluorides.
Outside oxides/fluorides, where no U is ever applied (n_elements=4):
  contains a U-list element        21.9%  [21.4, 22.4]  n=29,482
  contains only Sc/Ti/Cu/Zn 3d      7.2%  [ 6.8,  7.6]  n=14,649
  no 3d at all                      7.4%  [ 7.2,  7.6]  n=54,795
Same pattern at n_elements=3 (28.4% / 10.8% / 9.6%).

Interpretation: the driver is open-shell, magnetic d chemistry, not the
correction scheme and not "transition metal" as a category. Sc/Ti (d0) and
Cu/Zn (d10) behave like main-group elements. This survives the reviewer's
first objection.

## Robustness

- Threshold curves 0-50 meV/atom: the ordering (oxide > pnictide >
  boride/carbide/silicide > chalcogenide ~ halide ~ hydride > intermetallic)
  holds at every threshold. Not an artefact of the 0 cutoff.
- Both reference hulls (MP and MP+OQMD) give the same ordering; OQMD rates
  run 2-4pp higher throughout.
- Newer-PAW-potential elements (Li/Na/Mg/Ge/Ga) flip LESS (11.0% vs 13.5%),
  so the potential change does not manufacture the effect.

## Still outstanding before drafting
1. Selection effect: 139,541 of ~381,000 recomputed. Diff MaterialIds against
   the GNoME summary CSV; report recomputation coverage per family.
2. Known-materials baseline (Q2) via Materials Project API.
3. Budget model (Q3).
4. Re-run the scoop check immediately before submission.

---

# Phase 3a: selection effect (run 2026-10-09, same session)

Population file: stable_materials_summary.csv, 554,054 unique MaterialIds.
r2SCAN file: 139,541 unique MaterialIds. 117,043 match the summary (21.1% coverage);
22,498 r2SCAN IDs are NOT in the current summary (flip 11.4% vs 13.4% for matched),
presumably entries recomputed then displaced off the released stable set.

## Finding 1 (near-total gap): quinary materials are essentially unvalidated

CORRECTED 2026-10-09 after a direct check. The earlier draft of this section said
"no quinary or senary material was recomputed at all". That overstated it.
Exact coverage by element count:

  n_elements   population   recomputed   coverage
      2            1,741          224      12.9%
      3           63,441       20,187      31.8%
      4          306,688       96,403      31.4%
      5          177,392          229       0.13%
      6            4,792            0       0.00%

So quinaries are not literally zero: 229 of 177,392 were recomputed, a 240-fold
drop in coverage versus quaternaries. Senaries are zero. Individual quinary cells
with n_pop >= 500 all show 0 recomputed:
  5-element intermetallic + Ln          n_pop 39,940   recomputed 0
  5-element intermetallic + Ln + 3d     n_pop 23,866   recomputed 0
  5-element boride/carbide/silicide+Ln  n_pop 21,517   recomputed 0
  6-element cells                       all 0

The r2SCAN subset is restricted to <= 4 elements. GNoME's headline novelty claim is
discovery in "combinatorial spaces of 4+ unique elements" (Fig 2a) and emergent
generalization to 5+ elements. That chemistry has ZERO functional validation.
This is reportable as a finding, not only as a limitation.

## Finding 2: coverage declines toward the hull, where disagreement concentrates

  PBE decomposition energy   coverage
  (-inf, -0.2]                 18.3%
  (-0.2, -0.1]                 25.8%
  (-0.1, -0.05]                26.1%
  (-0.05, -0.02]               23.6%
  (-0.02, -0.01]               20.1%
  (-0.01, 0.0]                 16.5%
  (0.0, inf]                    0.3%

Supplementary Note 2 states the PBE/r2SCAN discrepancy concentrates near the hull.
The recomputed subset under-samples exactly that region, so the naive flip rate is
biased DOWNWARD. Report 13.0% as a lower bound for the <=4-element population.

## Finding 3: reweighting

Naive flip rate (recomputed set)        13.04%
Reweighted to the <=4-element population 11.96%  (63 cells, 384,242 materials)
Censoring bounds over the full 554,054: 3.3% (none of the unrecomputed flip) to
78.1% (all flip). Bounds are uninformative; reweighting is the defensible estimate.

## Consequence for the paper

The claim must be scoped: "among GNoME materials with <= 4 elements that DeepMind
recomputed at r2SCAN". The chemistry-resolved spread (26.7pp) is unaffected, because
coverage within the <=4-element cells is uniform at 41-46% with no strong gradient.
The GO verdict stands.

Outstanding: known-materials baseline (Q2), budget model (Q3), and an explanation for
the 22,498 unmatched IDs.

# Phase 5: verification-budget model (Q3)

Cell flip rates fitted on a random half of the recomputed set, catch evaluated on
the held-out half (n=69,771, 9,078 wrong verdicts, 13.01%). Ranking by rates fitted
on the same data would be circular.

  budget   chemistry-ranked   random   lift
     5%          11.5%          4.9%   +6.6pp
    10%          22.0%          9.9%  +12.1pp
    20%          38.1%         19.9%  +18.2pp
    30%          51.6%         29.9%  +21.7pp
    50%          71.6%         49.8%  +21.9pp

At a 20% recomputation budget, chemistry-aware allocation catches 1.9x as many
wrong stability labels as uniform sampling. This is the actionable contribution.

# The 22,498 unmatched IDs: benign

Same composition profile as the matched set (intermetallic / pnictide /
boride-carbide-silicide in the same order; n_elements distribution proportional),
flip 11.4% vs 13.4%. Most likely a release-version mismatch between the r2SCAN file
and the current summary, not a distinct population. Including or excluding them does
not change any conclusion. Still worth one line in the paper.

# Phase 4: known-materials baseline (Q2) — run 2026-10-09

Source: Materials Project thermo endpoint. NOTE the thermo_type value is "r2SCAN"
(case-sensitive); "R2SCAN" silently returns total_doc=0.
Pulled 37,123 r2SCAN docs and 152,426 GGA_GGA+U docs; 35,176 material_ids carry both;
22,756 of those are GGA-stable (energy_above_hull <= 1e-6) and so can flip.

Known-materials flip rate, aggregate: 14.88%
GNoME flip rate, aggregate:           13.04%

## Taken at face value, GNoME looks BETTER than known materials. Stratified, it reverses.

20 cells match between the two sets. GNoME is higher in 14, lower in 6, CIs disjoint in 9.

  cell                              GNoME    known    delta
  4el oxide, 3d, no Ln              30.7%    13.3%   +17.3pp  (disjoint)
  3el chalcogenide, 3d, no Ln       30.3%    17.8%   +12.5pp  (disjoint)
  3el boride/carbide/silicide, 3d   29.6%    18.3%   +11.3pp  (disjoint)
  3el halide, no Ln/3d              17.3%     9.2%    +8.2pp  (disjoint)
  4el oxide, no Ln/3d               18.9%    10.8%    +8.2pp  (disjoint)
  3el oxide, no Ln/3d               21.7%    16.7%    +5.0pp  (disjoint)
  3el pnictide, no Ln/3d            11.1%     6.9%    +4.1pp  (disjoint)
  ...
  3el intermetallic, Ln, no 3d       6.4%    13.0%    -6.6pp  (disjoint)
  3el intermetallic, no Ln/3d        7.9%    15.7%    -7.7pp  (disjoint)

Reading: AI-discovered materials are markedly LESS reliable than known materials in
oxide and 3d-containing chemistry, and MORE reliable in simple intermetallics. The
aggregate comparison hides the reversal because GNoME's composition mix is dominated
by intermetallics, where it does well. This is the strongest result in the project:
Merchant et al. compared 84% against MP's 90% in aggregate and concluded comparability.
Stratified, comparability fails in both directions.

## Limitations to state
- MP "known" includes theoretical entries, not only experimentally realised ones.
- MP prioritised r2SCAN for near-hull and <=20-site materials, so the known baseline
  is itself a selected subset; same caveat as the GNoME side.
- 4-element known cells are small (n = 108-297). Treat those deltas as indicative.
- Both sides use each functional's own hull; neither is ground truth.

# Phase 6: referee rounds, hostile pass and reproduction check — 2026-10-10

Three referee reports were generated and answered (REFEREE_REPORT_3.md in the paper
build is the hostile pass). Analyses added in response, all run from results/:
- phase7/phase8/phase10: prototype clustering (8,429 clusters, ICC 0.138, design effect
  3.14), cluster bootstrap CIs everywhere, BH over elements, budget over 20 cluster-wise
  splits, sandwich SEs on the logit (2.77x inflation on open-shell d).
- phase9: propensity/IPW for recomputation selection. Spread 22.0 -> 21.9pp. Null.
- phase11: pre-registered gate 6 (hull-distance strata) had NOT been run. Run late;
  gradient survives in every PBE band. Directional shift Delta added (oxide+3d +11.3,
  pnictide+3d +9.6, intermetallic -12.4 meV/atom). Fe/Co excess by anion class
  reverses the Section 4.3 framing: +5pp in oxides (U applied) vs +34 to +41 outside
  (no U). 25 non-PBE-stable rows excluded.
- phase2 "All" reference promoted to partial hull-completeness probe: oxide/intermetallic
  ratio 2.1 -> 1.4, ordering preserved.

Reproduction: run_all.py --check compares the 12 abstract quantities against results/.
12/12 matched on 2026-10-10. The check found that phase10_robust.py's last logged run
had crashed after the regression (KeyError '_g'); the script was fixed and re-run, and
the alternative-clustering intervals quoted in Section 3.3 reproduced exactly.

Open at this point: MP-side r2SCAN hull-completeness probe (q2b_completeness.py),
repository URL, endorser, pre-submission scoop check.

# Phase 6b: MP-side hull completeness (q2b_completeness.py) — run 2026-10-10
Median r2SCAN completeness by anion: oxide 37.0%, bcs 42.0%, pnictide 53.3%, halide 58.2%,
chalcogenide 58.5%, intermetallic 64.1%, hydride 67.3%. 3d excess by completeness tercile:
+8.1 / +9.9 / +5.3pp. Not a completeness artefact. Hull-distance matching vs MP impossible
(energy_above_hull floored at 0). Corrected the paper's claim about MP hull construction:
the r2SCAN thermo type is NOT restricted to fully-recomputed systems; the mixed scheme is.

# Phase 6c: CI caught a provenance gap — 2026-10-10

The first clean-runner reproduction (GitHub Actions, run 38042500440) came back green
with phase11_hostile reporting FAILED. Two separate defects:

1. phase11_hostile.py crashed at the rank-correlation line: pandas .corr(method="spearman")
   imports scipy, which is deliberately not a dependency. Both series were already ranks,
   so the argument was redundant; removed. The crash was silent in earlier sessions because
   sections A-C print before it.
2. Everything after that line had therefore never executed. Section D (within-anion
   Fe/Co/Ni excess) is the source of Table tab:feexcess in Section 4.3 — the table the
   round-3 referee pass called the paper's most important reframing. The numbers in the
   paper were correct but were produced ad hoc and no committed script regenerated them.
   Section D rewritten to compute Fe, Co and Ni excess over the no-that-element baseline
   within each anion class, writing results/within_anion_excess.csv. All 21 rows reproduce
   the paper exactly (Fe oxide +5.0, bcs +41.0; Co oxide +4.6, bcs +34.0; Ni oxide +10.3,
   intermetallic +2.0).
3. run_all.py now exits non-zero when any phase fails, and the workflow deletes the
   committed results/ before running, so a green badge cannot rest on stale files.
   Eight Section 4.3 quantities added to the headline check (now 20 checks, not 12).

Lesson recorded because it generalises: results committed to the repository make a
reproduction check self-confirming unless they are deleted first.

The first run under the new guard went red, correctly: deleting results/ also removed
known_baseline.csv, which the Materials Project phase produces and which the runner
cannot regenerate without an API key, so figures.py crashed. All 20 checked quantities
reproduced from the raw CSVs in that same run. The two MP-derived files are now exempt
from the delete and any check reading them is labelled "[stored, not regenerated here]".

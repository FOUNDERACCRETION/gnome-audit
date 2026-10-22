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

# Phase 7: independent adversarial review and major revision — 2026-10-10

A cold referee that had not seen the drafting was asked to review Draft 7 to reject. It
returned major-revision-or-reject. Everything it raised was checked against the data
rather than accepted; the outcome was five real defects, two of which changed claims.

CONFIRMED AND FIXED
1. BH correction was broken: cummin applied over an ascending p-value sort propagated
   Fe's p~1e-85 to every later element. Hf had raw p=0.95 and was flagged significant.
   The abstract claimed all 83 contrasts survived. True count: 51 of 83.
2. The budget model benchmarked chemistry-ranking only against uniform random. Ranking
   by PBE hull distance alone — free, no fitting — beats it at every budget (40.4% vs
   38.1% at a 20% budget). Combined reaches 46.1%. The paper's operational
   recommendation was wrong and is now "hull distance first, composition on top,
   worth ~6 points".
3. Sections 3.6 and 4.10 stated OPPOSITE pre-registered lanthanide predictions, and both
   were backwards relative to the committed spec, which registered "at or below the
   global rate, refuted if above". Result 11.8% vs 13.0% CONFIRMS it. Fixed.
4. "Borides, pnictides and intermetallics, where PBE runs with no correction at all" is
   false: MP2020 corrects O S F Cl Br I N H Se Si Sb Te. Only the +U is cation-specific.
   Reframed around that, which is cleaner: anion corrections are per-anion constants and
   cancel in a within-anion excess, so the design already isolated U. Added the decisive
   contrast: within halides, F-containing (U applied) Fe excess +8.8 vs +16.3 for
   Cl/Br/I only; Mn +2.4 vs +21.5.
5. The claim that MP exposes no sub-hull margin for stable entries was wrong (pymatgen
   get_decomp_and_phase_separation_energy recovers it). Stated as not-yet-done rather
   than impossible.

REFEREE POINTS CHECKED AND REBUTTED
6. "Elemental reference offsets could produce the whole gradient." They cannot: a
   decomposition energy is evaluated at fixed composition, so an additive per-element
   shift moves compound and competitors identically and cancels. Reference errors move
   formation energies only. The empirical version of the test is also collinear by
   construction, because the cells are defined by element presence — noted in the paper
   so no one repeats it. The discriminating test is context dependence: Fe contributes
   -31 meV/atom per unit atomic fraction in oxides and +195 in halides, spread 239, where
   additivity requires 0. Cr 249, Mn 273. For Fe/Cr/Mn the oxide context is lowest, an
   independent recovery of the U-suppression result.
7. "Winner's curse." Gradient is +20.0 vs +22.3 across generation pipelines and +21.6 vs
   +19.4 across training-set membership. Rules out a differential curse; a uniform one
   remains possible and is now stated as a live alternative for the MP comparison only.

OTHER CORRECTIONS
8. Validation compared mismatched populations. Quaternary survival 87.35% vs reported
   86.8% (0.55pp); binaries+ternaries 85.34% vs 84%. The "within 0.2pp" claim was a
   coincidence between different populations and is withdrawn.
9. Design effect: equal-size formula gives 3.15, unequal-size gives 46.3, because five
   clusters hold 7.7% of the data. Neither is used for inference; both are now reported,
   with the empirical bootstrap (1.82) and sandwich (2.77) landing between. Dropping the
   5/20/50 largest clusters moves the headline from +26.7 to +26.8/+27.3/+28.2.
10. Band gap: the old binary split concluded the gap dependence "largely disappears"
    once 3d is controlled. False. In bins, the 3d excess is +10.5/+9.6/+10.4/+11.8 and
    the gap effect persists at fixed 3d content (8.5% to 18.9% among no-3d). Two
    independent axes. Also: the gap column is in the summary file, not the r2scan file,
    and contains infinities that were never screened.
11. Logistic model was main-effects-only. The lanthanide OR of 1.02 averages -9.9pp in
    halides against +7.0pp in pnictides; interactions improve fit by chi2=879 on 14 df.
12. Shift table rebuilt on the primary cell definition with quartiles and 1% shares,
    fixing cell-count discrepancies of 2-10 materials against other tables.
13. Abstract cut 545 -> 378 words; draft header comments stripped; repeated disclaimers
    and the "nor is it X" section cadence removed; "All"-reference inference explicitly
    marked as non-diagnostic; benchmark paragraph narrowed to what the data support.

run_all.py now checks 23 quantities including two claim-level checks: that the budget
ordering holds at every budget, and that the Fe excess is smallest in oxides.

# Phase 8: second cold error-hunt — 2026-10-11

A second independent agent, told to hunt defects rather than review, returned 74 on
Draft 10. Each was checked against the data before acting. Five were rated FATAL and
four of those were real.

FATAL, fixed
1. "13.0% should be read as a lower bound on the population rate" contradicted Section
   4.13: IPW moves the estimate DOWN, 13.37% -> 12.71%. The hull-distance component does
   push up but is outweighed by composition and cell size. No direction is now claimed.
2. The regressor labelled "open-shell d (U list)" IS the U list, which contains Mo and W
   (4d/5d) and excludes Cu (d9 as Cu2+). The paper then argued the operative variable was
   magnetic open-shell d character rather than U-list membership, which that regressor
   cannot test. Refit three indicators (phase21): U list logLik -50471.5, 3d subset
   -50328.5 (best), open-shell 3d incl. Cu -50589.5 (worst). Cu flips at 12.6%, BELOW the
   global 13.0%, despite being open-shell and magnetic; Mo and W flip at 20.2% and 18.2%.
   The magnetic reading is refuted by our own data and has been replaced: the predictor is
   U-list membership, i.e. the property that motivated the correction list.
3. Delta was described as a functional difference. It is (uncorrected r2SCAN) minus
   (GGA+U + MP2020), so it contains the correction removal. Defined explicitly now, and
   the element x anion table is no longer called an "independent estimate" of the U
   suppression -- it is the same measurement reparameterised.
4. The pre-registered lanthanide test compared a group that is 59.4% of the dataset
   against a global rate it dominates. Against non-lanthanides the gap is 11.8% vs 14.8%,
   and the regression-adjusted effect is null (OR 1.02 [0.93, 1.12]). Reported as
   confirmed-as-stated but uninformative as posed.
5. NOT real as stated, but exposed our own error: the band-gap section claimed 26,995
   materials carry an infinite gap. The true figure is 2,177 infinite; the attrition is
   dominated by the 22,498 IDs absent from the summary file. Rewritten.

SERIOUS, fixed: the "All" ordering is not preserved in detail (chalcogenides fall below
intermetallics; oxides no longer rank first); intermetallics gain most only on the ratio
scale; hydrides not intermetallics are nominally the most complete MP class; Ni does show
a small intermetallic excess (+2.0 [+0.3,+3.8]); the halide interval excludes unity;
anion corrections do not cancel exactly and the pnictide and boride classes are mixed with
respect to them; "PBE" now defined once as PBE+U+MP2020; Table 13 gained n, cluster
intervals and a sign test (14/20 worse, p=0.115, 11 of 20 intervals exclude zero); Fe
excess gained an odds-ratio column since oxides have the highest baseline.

One apparent contradiction worth recording because it will be raised again: fluorides
flip at 24.8% against 7.2% for Cl/Br/I, and F is the 5th-highest element overall, which
looks like it refutes U-suppression. It does not -- the class baseline and the
within-class TM excess are different quantities. The excess is +6.3 among fluorides
(where U applies) and +7.8 among the other halides. Both are now in Section 4.3.

Remaining open: the budget comparison still lacks a paired per-split difference; the
Matbench test-set composition is still uncomputed; the "All" column is still undocumented
and unasked; references still not author-verified.

# Phase 9: third error hunt, reproducibility hardening — 2026-10-11 (overnight)

Third independent hunt on Draft 11 returned 41 defects. Each checked before acting.

FATAL, fixed
- The hull-completeness bias on the known-materials comparison was stated with the sign
  inverted: if the MP oxide flip rate is understated then GNoME-minus-known is OVERstated,
  so +17.3 is an upper bound, not a lower bound. Introduced during the round-4 rewrite.

SERIOUS, fixed
- Appendix B claimed "all 83 elements" but the source table contained only 59, and the
  multicols environment silently overflowed the page. Rebuilt from
  per_element_clustered.csv as a longtable: all 83 now render, with an asterisk marking
  the 51 that survive BH. Copper, the element the Section 4.3 argument turns on, was one
  of those missing from the printed table.
- "Oxides, pnictides and borides remain highest and intermetallics remain lowest" under
  the All reference is false: chalcogenides (21.6%) fall below intermetallics (22.3%) and
  oxides drop to third. Narrowed to the top-versus-intermetallic contrast.
- Section 4.2 still carried the pre-round-4 reading (open-shell character is the operative
  variable) that Section 4.3 now refutes, and misdescribed Ni's pattern. Rewritten.
- The U-suppression claim was stated for the U list as a whole; it holds for Fe and Co.
  Ni's largest excess is in oxides (+10.3), i.e. the opposite direction. Now stated as a
  Fe/Co result with Ni named as the exception.
- "The ordering is unchanged on either scale" was contradicted by the odds column: on
  odds, intermetallics (8.38) overtake borides (7.88). Narrowed to the bottom two.
- "Every cell on which GNoME is significantly worse contains oxygen, a 3d metal or a
  halide" is false for two of eight (ternary pnictides and chalcogenides, the two
  smallest deficits).
- Two adjacent contradictory sentences in Section 4.6 ("The cells are not differently
  positioned" immediately followed by "The cells do differ somewhat in position").
- Copper was called "the decisive case" on an assumed Cu(2+) oxidation state in a dataset
  that is >60% metallic, where most Cu is likely d10. Downgraded to suggestive.
- Lanthanide-pooled cells in Tables 6 and 7 and Section 4.5 were never flagged as pooled,
  making the headline contrast appear to shrink from 26.7 to 22.0 to 20.0 across sections.
  Captions now state the pooling.
- Nine places still said "open-shell d" after the label was changed to "MP U-list".
- LR test df was 14; with five anion dummies it is 10.

CITATION ERROR FOUND AND REMOVED
- rare2026 (Tao, Ridwan, Ke, Zhu, arXiv:2608.25270) was cited for the claim that GGA
  treatment of the Gd 4f7 shell renders the MP hull unreliable. Verified at source: that
  paper does not mention Gd, the 4f shell, GGA, or MP hull reliability, and Gd is not
  among the elements it studies. Claim and citation deleted. warford2026 verified at
  source and its author initials completed. liu2024 verified as a non-peer-reviewed
  commentary and now labelled as such at both points of use.

REPRODUCIBILITY
- scripts/check_paper.py: parses every number out of the LaTeX tables and diffs it
  against results/. 232 values checked. It found a real transcription error on its first
  run (MP completeness intermetallic n printed 7,953 against 7,954 in the data), which
  also explained a 22,755-vs-22,756 total that had been papered over with a footnote
  about an unclassifiable composition. Now wired into run_all.py and the CI workflow, so
  table-vs-analysis agreement is machine-checked on every push rather than trusted.
- Budget comparison gained paired per-split differences, which is the correct uncertainty
  for policies evaluated on shared splits: hull beats composition by +2.24 pp (paired sd
  1.09, 20/20 splits) and combined beats hull by +5.77 pp (paired sd 0.95, 20/20).

# Phase 10: fourth external review — stale figures found — 2026-10-11

A fourth reviewer (7/10) found what no automated check could: the figures compiled into
the PDF were stale copies in the build directory, not the regenerated ones in results/.
Three were wrong:
- Fig 3 (band gap) predated the Section 4.8 rebinning: it showed a "metal (0)" bin with
  n=79,081 that no longer exists in the analysis.
- Fig 4 (logit) plotted naive intervals and raw variable names (has_U, is_bcs).
- Fig 6 (budget) showed only two of the four policies and was labelled "wrong stability
  labels caught", contradicting the paper's disagreement-not-error framing.
All three regenerated from current results (scripts/figures5.py for 3 and 4,
figures4.py for 6) and copied into the build. check_paper.py cannot catch this class of
error because it reads tables, not images; a figure-provenance check would need the
plotting scripts to stamp their inputs, which is left as future work.

Also fixed from that review:
- Tie-breaking in the budget evaluation was a stable sort on input order. Where the
  budget cuts through a tied cell this biases the realised catch. Now a seeded random
  permutation. Numbers shift slightly (composition 38.1 -> 37.6, combined 46.1 -> 45.8
  at a 20% budget); all orderings and all 20/20 paired wins survive. Paired differences
  are now hull-composition +2.75 (sd 1.43) and combined-hull +5.42 (sd 1.39).
- Causal U language removed throughout: the comparison is across composition classes
  that differ in more than whether U is applied, so the paper now says the smaller oxide
  excesses are consistent with partial suppression, and explicitly declines the claim
  that applying U causes the gradient to shrink.
- "Upper and lower bounds" on the chemistry effect replaced with reference-set
  sensitivity. Two observed reference choices do not bracket the effect; a more complete
  hull raises individual decomposition energies but need not move class rates equally.
  The 2.1-to-1.3 ratio range is now presented as a result rather than a caveat.
- "Excess attributable to Fe" -> "associated with"; "same chemistry" -> "same broad
  composition class".
- Hull stratification no longer claims to "remove" proximity differences; binning leaves
  residual variation and a continuous model is named as the stronger untried test.
- The F-statistic and the LR chi-squared are now labelled as assuming independence and
  read descriptively, consistent with the paper's own clustering correction.
- The 20 matched MP cells are labelled exploratory with no multiplicity correction.
- The budget is called an entry-count budget, not compute.
- Section 4.13 said the aggregate inversion was "consistent with" the source's 84-vs-90
  comparison. It is opposite in direction; now stated as such.
- Section 4.12's "the expectation fails" clarified: the frozen-f concern is not borne
  out, our registered expectation is.
- Citation-use fixes: the introduction's blanket "every one of those verdicts is PBE" is
  false (Alexandria uses PBEsol geometries and SCAN energies) and is now scoped to GNoME;
  Kirklin's 0.081 qualified as the all-elements fit (0.096 distributed) and the 0.082
  experiment-vs-experiment figure identified as 75 overlapping intermetallics rather than
  a universal noise floor; PseudoDojo cited as the general statement with an explicit
  note that GNoME's actual lanthanide PAW datasets are unknown to us; Zeni completed to
  Nature 639, 624-632; Reents given author initials; the MP mixing documentation now
  cites the specific page rather than the root domain.

# Phase 11: fifth review — population accounting — 2026-10-11

A fifth reviewer (major revision) returned ~120 recommendations. Triaged: about eleven
were real defects, the rest largely hedging preferences, several already done. Fixed the
eleven; declined most of the rephrasing, which applied wholesale would leave the paper
making no claims at all.

The important one: the band-gap population fix reported as done in Phase 9 HAD NOT
LANDED. It was written, the edit batch aborted on an unrelated line, and the re-applied
batch omitted it. The paper still said 26,995 materials carry an infinite recorded gap.
phase22_populations.py now produces the mutually exclusive accounting, which the paper
carries as a table: 22,498 unmatched / 2,320 empty field / 2,177 infinite sentinel /
112 exactly zero / 112,434 finite positive = 139,541. Only the 2,177 are sentinels and
none of the groups is evidence of insulating character. The 112 exact-zero figure is
verified. Lesson: check_paper.py validates tables against data but nothing validates
prose against data, and a reported-fixed item was wrong for two days as a result.

Other real fixes:
- Section 4.3's heading still read "The correction scheme suppresses the disagreement
  where applied", contradicting the cautionary paragraph added in Phase 10, as did
  "Everything below therefore isolates U", "the empirical correction is doing part of
  what the meta-GGA does" and "neither surface is corrected". All replaced. Paper-wide
  formulation is now: elevated disagreement persists outside the U-treated regime; the
  smaller Fe/Co excesses in oxides are consistent with partial suppression but do not
  identify a causal effect.
- Section 5 said larger-celled materials flip less often; Section 4.10 says rates rise
  from 11.3% to 16.3% with cell size. Contradiction removed and the per-feature
  decomposition of the IPW direction withdrawn, since the model as fitted does not
  support it.
- Nickel's largest excess is chalcogenides (+14.2), not oxides (+10.3).
- Censoring bounds were computed treating all 139,541 recomputed entries as observed
  members of a release 22,498 of them do not match. Recomputed on the eligible
  population: N=365,504 PBE-stable <=4-element release entries, n=116,789 recomputed,
  F=15,608 flips, giving 4.3% to 72.3%.
- "The four lowest" listed five elements; Figure 1's dashed line is the all-element rate
  not the quaternary marginal; Kuner completed to npj Comput. Mater. 11, 209.
- Added Table 1, a population-flow table, so every varying n in the paper traces to a
  stated exclusion.
- Abstract: "correcting throughout" made specific; 239 meV/atom identified as a fitted
  coefficient per unit atomic fraction; "All" no longer called a more complete hull;
  scope narrowed to the primary cell analysis; "unvalidated" replaced with "extremely
  sparsely recomputed".
- Introduction: the benchmark sentence no longer asserts what benchmarks assume, and a
  paragraph now separates formation-energy accuracy from verdict disagreement.

Declined: the wholesale replacement of every "X does not generate the effect" with
"persists under the tested adjustment for X"; the full provenance list in the reviewer's
6.2 beyond access dates and query parameters; deleting the F-statistic already labelled
descriptive. The reviewer's own section 13 stopping rule is adopted: remaining work is
the GNoME authors' email about the "All" column and an independent human check of the
core scripts.

# Phase 12: systematic consistency pass — 2026-10-11

Sixth review (8/10, major revision). Its diagnosis was structural rather than a list:
the body still carried conclusions the revised abstract had dropped. The cause is
patch-based editing — passages reviewers flagged got fixed, passages nobody flagged were
never touched. Section 1.2, the contributions list, had not been edited in sixteen
drafts and still read "controlled tests", "is not produced by", "partially suppresses
rather than creates it". It also still contained a grammatical break introduced by the
previous patch round, which is evidence for the same point: patching had reached the
rate where it introduced about as many defects as it removed.

Response was one systematic pass rather than another patch batch: the abstract was taken
as the authoritative statement of the paper's claims and the contributions list,
discussion, limitations and conclusion were rewritten as a set to agree with it. The
conclusion was rewritten entirely.

Substantive changes beyond wording:
- The element-versus-global bootstrap resampled the element rate alone against a fixed
  global rate, although the comparator is estimated from the same data. Now resampled
  jointly: one global cluster draw per replicate, both rates recomputed from it, standard
  error taken from the difference. The BH count moves from 51 to 53 of 83. Appendix B
  asterisks regenerated; run_all check updated.
- Hull-coverage claims narrowed: entry count is a proxy and does not determine the
  magnitude or ordering of the resulting bias, so the claims that oxide rates are "most
  understated" and that completing the hull "would widen" the contrast are withdrawn.
  What survives is that the 3d contrast persists within the tested coverage strata.
- The winner's-curse exclusion is withdrawn. A common energy displacement need not
  produce equal changes in flip probability across classes, so "a common shift cancels in
  a difference" was wrong.
- The flip outcome is now defined on two populations explicitly: confirmed
  PBE-stable-to-r2SCAN-unstable for the 117,043 matched materials, and positive r2SCAN
  decomposition energy for the 22,498 unmatched where PBE stability is inferred rather
  than observed.
- The IPW target population (371,870, all <=4-element entries) is distinguished from the
  censoring-bound population (365,504, those also PBE-stable); they differ by 6,366.
- McFadden pseudo-R2 no longer described as variance explained.
- "A genuine open-shell 3d set" is now "an alternative 3d element-presence indicator",
  and the conclusion drawn is that this indicator fits worse, not that open-shell
  character is excluded as a mechanism.
- New Appendix D, implementation details: resampling unit and replicate counts, joint
  contrast resampling, regression indicator encoding including the Mo/W and Fe+Ti cases,
  sandwich estimator with residuals and scores defined separately, propensity weighting
  with its positivity problem stated, and the budget evaluation's eligible cluster count
  (7,822), filter-before-split order, generated-regressor caveat and analytic
  uniform-random row.

Both checkers green: 23 headline quantities, 232 table values.

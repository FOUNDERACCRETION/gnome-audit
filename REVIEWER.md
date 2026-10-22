# Independent review packet

This file exists because the one open item on this project cannot be closed from inside
it. Two separate checks are needed and they are not the same job.

**Check A — can you reproduce it?** Can someone with the inputs and this repository
regenerate the central tables and figures on their own machine, without using the stored
results?

**Check B — is the analysis right?** Running the scripts proves they run. It does not
prove the eligibility rules are the right rules, the encodings encode what the paper says
they encode, or the contrasts are resampled correctly. Check B is reading, not running.

Both have already caught things. Of the defects found so far in this work, **every single
one was found by reading code against the manuscript, not by running it.** Both automated
layers were green through all of them:

| Defect | How it was found |
|---|---|
| Benjamini–Hochberg `cummin` applied over an ascending p-sort, propagating Fe's p ≈ 1e-85 to every later element | reading the step-up code |
| A whole section of `phase11_hostile.py` had never executed, after an exception on a `.corr(method="spearman")` call | reading a log that nobody had read |
| Hydride sensitivity fitted with hydride as its own anion dummy, so its baseline was not the regression table it claimed to test | a reviewer noticing 1.103 ≠ 1.101 in print |
| IPW weights documented as untrimmed while the code floored the propensity at 1e-2 — which turned out to change the estimate from 3.05% to 12.71% | reading the weighting code |
| Bootstrap endpoint instability explained away as a property of the sampling distribution rather than as Monte Carlo error | a reviewer applying the definition of a percentile estimator |
| The matched-comparison intervals drawn from independent binomials on both sides, so the GNoME side carried none of the prototype clustering the paper claimed for its rate intervals. Correcting it widened them 1.28x and moved one cell across zero | a reviewer noticing that "binomial resample" and "cluster bootstrap" cannot both be true |

That is the calibration for how much weight to put on the green checks. They catch stale
outputs and transcription errors. They do not catch a wrong definition.

---

## Check A: reproduction

Follow `REPRODUCE.md`. The short version:

1. Python 3.10+, `pip install -r requirements.txt` (pandas, numpy, matplotlib; no scipy,
   no scikit-learn, no pymatgen).
2. Download the two GNoME CSVs into `data/`. Verify the SHA-256 sums printed in
   `REPRODUCE.md` §2 — not just the row counts.
3. **Delete `results/` first**, except `known_baseline.csv`, `mp_completeness.csv`,
   `mp_gga_raw.csv` and `mp_r2scan_raw.csv`, which come from the credentialed Materials
   Project pulls. Otherwise the checks can pass on files you did not generate. This is
   what the CI workflow does, and it is the single most important step.
4. `python run_all.py`. It regenerates `results/`, checks 54 headline quantities against
   the values printed in the paper, then runs `scripts/check_paper.py`, which parses every
   number out of the LaTeX tables and diffs 242 printed values against the result files.
   Per-phase wall times land in `results/runtimes.csv`.
5. Rebuild the PDF: `cd paper && pdflatex gnome_audit_paper.tex` three times.

A clean run ends with `ALL HEADLINE NUMBERS REPRODUCED` and `every parsed table value
matches the analysis output`. Anything else is a finding, including a phase that crashes
while the surviving numbers still match — the driver exits non-zero on that case
deliberately, because it has happened.

If you would rather not take the pipeline's word for anything, the most useful thing you
can do is **reimplement the headline composition-cell rates and the budget-policy
comparison in a short script of your own**, from the two CSVs, without reading ours. That
tests the logic rather than re-executing the implementation.

---

## Check B: is it right?

What follows is ordered by how much of the paper collapses if the answer is no.

### B1. The outcome and the eligible population — `scripts/phase2.py`, §3.3–3.4

- `flip` is `Decomposition Energy Per Atom MP > 0`. Is a strict inequality right, given
  that the field is a float and the hull is defined by equality?
- The paper defines two outcomes: *release-based disagreement* (all 139,541 rows, no
  reference to the PBE verdict) and a *confirmed transition* (the 117,018 rows matched to
  the summary file with a non-positive PBE energy). Check that every number in the paper
  uses the one the text claims, especially §4.13, where the Materials Project side is a
  confirmed transition and the GNoME side is release-based. That asymmetry is stated, and
  the comparison is also reported on the matched basis — verify both.
- 25 matched entries were not PBE-stable. The paper retains them in release-based rows and
  they are outside the confirmed-transition populations by construction. Verify rather
  than assume.
- 22,498 r²SCAN IDs have no summary row. Their PBE stability is *inferred* from release
  provenance. Is the paper's treatment of them consistent everywhere?

### B2. The clustering key and everything resting on it — `phase7_cluster.py`, `phase34_converged_intervals.py`, §3.5

Every interval in the paper is a cluster bootstrap over this key, so an error here moves
all of them.

- The key is `(space group number, anonymised stoichiometry)`. Compositions are used **as
  released and not reduced**, so the key depends on unit-cell convention. How much does
  that matter? A reduced-stoichiometry sensitivity check has not been run and is the
  obvious next one.
- The anonymiser maps a composition to sorted site counts (`Nd7Os1Pr3Si5` → `A1B3C5D7`).
  Read the regex. Does it handle every formula in the file, including single-atom
  stoichiometries written without a `1`?
- The resampling in `phase34` draws multinomial counts over clusters instead of drawing
  cluster indices. Satisfy yourself those are the same distribution, and that the
  per-cluster count and flip-sum aggregation is equivalent to resampling rows.
- ICC 0.138 and the design effects (3.15 equal-size, 46.3 unequal-size): check the
  formulas and which one the paper leans on where.

### B3. The group encodings — `phase6b_logit.py`, `phase10_robust.py`, §4.11 and Table 2

- `has_U` = any of V, Cr, Mn, Fe, Co, Ni, W, Mo. `has_3d_nonU` = one of Sc, Ti, Cu, Zn
  **and** no U-list element. The two are mutually exclusive by construction; the anion
  dummies pool intermetallic **and hydride** into the reference. Verify that Table 2's
  three groups partition the eligible set under those definitions — a Mo-only composition
  has no 3d metal but is in the U-list group.
- The hydride-pooling sensitivity (`phase29_reconcile.py`) must start from exactly the
  Table 10 design. It did not, once. Check that the "with hydrides" column reproduces the
  printed table to three decimals.
- MP2020 applies the Hubbard U **only in oxides and fluorides**. The whole §4.3 argument
  turns on that. Verify it against the Materials Project documentation, not against us.

### B4. The bootstrap contrasts — `phase8_referee.py`, `phase16_excess_ci.py`, `phase34`, `phase37_matched_cluster.py`

- The matched comparison (`phase37`) resamples the two sides differently and this is the
  part most recently wrong: prototypes with replacement on the GNoME side, drawn **once
  across all twenty cells per replicate** so that cells sharing a prototype stay dependent;
  material-level on the Materials Project side, because no prototype key exists for those
  entries. Check that the aggregation by (prototype, cell) is right, that the standardised
  contrast uses the same draw rather than a fresh one, and that the MP side's binomial draw
  really is the nonparametric bootstrap for a mean of binary outcomes. The superseded
  all-binomial intervals are still computed alongside, so the correction's size is visible —
  check that the paper quotes the clustered ones.

- The element-versus-global contrast resamples **one** cluster set per replicate and
  recomputes both arms from it. An earlier version resampled the element arm against a
  fixed global rate and gave 51 significant contrasts instead of 53. Check the current
  code does the joint version.
- The comparator includes the element under test, which attenuates the contrast for
  large-n elements. Is that the right comparator?
- Standard error → p-value is a two-sided normal approximation, then Benjamini–Hochberg
  over 83 tests. Check the step-up direction: `[::-1].cummin()[::-1]`, not `cummin()`.
- Replicate counts: 20,000 for cells and per-element, 100,000 for the within-anion
  excesses. `phase33_percentile_conv.py` measures the convergence. Check that the counts
  are high enough that the printed endpoints are stable to the precision they are printed
  to, and that the paper no longer claims any endpoint instability is irreducible.

### B5. The propensity model — `phase9_selection.py`, `phase26_diag.py`, `phase29_reconcile.py`, §4.9

This is the weakest part of the paper and the reviewer's stated priority.

- The model is fitted on 371,870 summary entries with ≤4 elements; the weighted outcome is
  computed on the 116,789 that are also PBE-stable and recomputed. Confirm the paper never
  conflates the fitting population with the outcome population.
- PBE decomposition energy enters linearly **and squared**, both standardised. Is that
  transformation stated identically in the model and in the diagnostics?
- The fitted propensity minimum over the eligible population is 7.1e-218. Is the IRLS
  converging, or is this separation and/or numerical failure? Inspect the fit, not just
  the output. The paper now reports this as estimated lack of overlap and explicitly does
  **not** claim population-level positivity failure — check that it keeps that line.
- Weights are floored at 1e-2. Untrimmed, Kish's effective sample size is **5** of 116,789
  and the estimate reads 3.05% instead of 12.71%. The paper reports 12.7% as a sensitivity
  result conditional on model and floor. Is that the right framing, and is the floor
  defensible?

### B6. The budget evaluation — `phase12_budget_baselines.py`, `phase24_budget_strong.py`, §4.14

- Splits are **by cluster**, 3,911 of 7,822 per split, eligibility applied before
  splitting. Verify no material's cluster straddles fit and evaluation.
- Ties within a cell are broken by a seeded random permutation, which matters because the
  budget cuts through cells.
- `combined_cf` cross-fits the cell rates in five folds **assigned to prototypes**, then
  mapped to materials. Verify the folds are cluster-wise and not material-wise.
- The uniform row is the analytic expectation `floor(bN)/N`, not a simulation.
- "Catch" = selected evaluation flips ÷ all evaluation flips. Check that definition against
  the code.
- The 20 splits overlap (each draws half the same pool), so the across-split spread is not
  a confidence interval. Check the paper says so.

### B7. The reference sets — §4.16, `phase20_audit.py`, `q2b_completeness.py`, `phase40_allref_cells.py`

- The `"All"` column is undocumented. Every conclusion depending on its identity is
  flagged conditional. Check that no claim leaks past that flag.
- MP r²SCAN "completeness" counts **every** entry in a chemical system and its subsystems,
  not just ground states, so it is an entry-coverage proxy. Check the paper never upgrades
  it to a measure of hull completeness.
- MP `energy_above_hull` is floored at zero (21.8% of 152,426 rows read zero), so sub-hull
  margins are not recoverable from it. Check nothing in the paper assumes otherwise.

### B8. Energies and references — §4.4

- A decomposition energy is invariant to a consistent composition-linear shift of all
  phases. Check the argument as stated, and check the paper distinguishes that from a
  change of elemental reference convention.
- Δ = E_dec(r²SCAN) − E_dec(PBE) contains the **removal of the MP2020 corrections**, not
  only a functional difference. Check every statement about Δ respects that.
- Table 4's model identifies row *spreads*, not levels (atomic fractions sum to one).
  Check the paper claims only the spread.

---

## What would make this a clean bill of health

Not "the code ran". Something closer to:

> I regenerated `results/` from the two CSVs on a clean machine after deleting the stored
> outputs, and the 54 headline quantities and 242 table values matched. Separately, I read
> the eligibility definitions, the clustering key, the group encodings, the bootstrap
> contrast construction, the propensity model and the budget splits against the
> manuscript's description of each, and I independently reimplemented [X] and got [Y].
> The disagreements I found are listed below.

A list of disagreements is a success, not a failure. The useful outcome of this review is
a finding, and the five defects above are the prior for how likely one is.

## Contact

Karthikeya — founder@accretionlabs.net. Please send findings as issues on the repository
if it is public by the time you read this, or by email if not. Everything you find will be
credited in the paper's acknowledgements unless you ask otherwise, and anything that
changes a number will be recorded in the response document in `review/`.

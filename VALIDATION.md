# Validation: what the checks assert, and what they cannot

`run_all.py` regenerates `results/` from the two GNoME CSV files and then applies three
layers of checking. `run_all.py --check` runs the checks against the stored results without
recomputing them. This file records what each layer asserts and why each check is
classified as it is; the manuscript states only the summary.

## Layer 1 — 84 headline quantities (`run_all.py`)

Compares 84 quantities against the values printed in the paper: the cell rates, the
regression odds ratio and its interval, the average marginal effect of U-list membership
with its bootstrap interval, the budget-policy ordering at every budget, the paired margins
of the combined policy over a nonparametric hull-distance model with and without
cross-fitting, the within-anion excess table, the analytical population counts, the
hydride-pooling sensitivity, the two headline composition cells under both reference-set
definitions, and the per-element multiple-comparison count with its stability across five
seeds.

Each check is labelled by kind, because the kinds carry different authority.

### `[impl]` — 6 implementation properties any legitimate run must satisfy

1. The baseline column of the hydride-pooling sensitivity reproduces the main regression
   table. A mismatch here means the sensitivity analysis is being run against a different
   model from the one it claims to perturb, which makes it meaningless rather than wrong.
2. The untrimmed maximum propensity weight equals the reciprocal of the smallest fitted
   propensity among weighted entries.
3. `phase38_resampling_equivalence.py`: the cluster-aggregated resampling path returns
   bit-identical cell rates to an explicit row-gather resample under the same cluster draw.
4. `phase38_resampling_equivalence.py`: a drawn prototype's multiplicity applies to every
   one of its members, in every cell in which it appears. This is the part that would fail
   silently if the aggregation were wrong, since 1,695 of these prototypes span more than
   one cell.

5. `phase43_context_identify.py`: the context-dependence design is full rank (113 of 113),
   so the reported entries are a unique least-squares fit and not an artefact of how the
   solver resolves degeneracy. The projection of each of the 75 within-element contrasts
   onto the numerical null space is exactly zero, the null space being empty.
6. `phase43_context_identify.py`: those contrasts are invariant under two invertible
   reparameterisations of the design, column scaling and column permutation, to
   6.5e-12 meV/atom. Note that dropping a column is *not* an invertible reparameterisation
   and an earlier version of this check wrongly used one; so did the first attempt at the
   scaling map-back, which must be `beta = D c`, not `c / d`. Both failures were in the
   check, not the model.

Checks 3 and 4 exist because the resampling aggregates to the cluster level before the
replicate loop. That optimisation is what makes 100,000 replicates affordable, and it is
the kind of optimisation that can be subtly wrong while still producing plausible numbers.

### `[result]` — 10 asserted findings

1. The budget ordering holds at every budget.
2. Composition adds after hull distance is modelled nonparametrically.
3. The combined policy's margin survives cross-fitting.
4. The iron excess is smallest in oxides of all anion classes.
5. The two headline cell intervals remain disjoint.
6. The oxide and boride iron intervals do not overlap.
7. The headline cell spread still excludes zero under the "All" reference set.
9. Four of the five elements reject equal interaction coefficients under
   Benjamini-Hochberg, and cobalt does not. Note what this replaced: an earlier version
   asserted that every context-spread bootstrap interval excludes zero. That is not a test
   of the additive null, because a spread is a data-selected maximum minus minimum, and the
   two disagree in exactly the predicted direction - cobalt's interval excludes zero while
   its joint test does not reject.
10. The U-list association survives continuous hull-distance adjustment.

These are findings, not invariants. If one fails on a legitimate re-run it is to be
re-examined and reported, not forced back into agreement.

### Unlabelled — 68 reference comparisons

Comparisons against a recorded reference run: staleness checks, expected to be updated
whenever a number legitimately changes. A passing reference check means a number has not
drifted from the run that produced it. It is not evidence that the statistical choice
behind the number was the right one.

One of these is easily mistaken for an invariant and is deliberately **not** labelled
`[impl]`: that the cluster-adjusted matched intervals are wider than the independence-based
ones in all twenty cells. That is a property of the dependence in these data, not a
requirement on any legitimate run — a cluster bootstrap can give a narrower interval than
an independence-based one for some statistics and dependence patterns.

Nothing in the set asserts that an interval excludes zero, and nothing asserts a rate of
Monte Carlo convergence. Both are results to be reported, not invariants for software to
enforce.

## Layer 2 — printed table values (`scripts/check_paper.py`)

Parses every number out of the LaTeX tables and diffs 655 printed values against the
result files, including the endpoints of every interval in the Fe within-anion excess
table, the matched-cell table, the two-reference headline cell table, the continuous
hull-distance table and the context-dependence table.

One trap worth recording: `\phantom{0}` alignment padding must be stripped before a cell
is parsed, or its zeros are read as data, and interval endpoints printed as integers need a
tolerance of half a unit, not a full one. A full unit silently accepts an off-by-one, which
a regression test caught. Neither a
transcription error nor a regenerated bootstrap can leave one of those tables stale.

Coverage is not automatic for intervals quoted in **prose** rather than in a table. Those
are enumerated one by one in `PROSE_INTERVALS` at the top of the script — the Ni
intermetallic excess is there because the prose once quoted an upper endpoint of 3.8
against the table's and the CSV's 3.7, and neither Layer 2 nor Layer 3 would have caught
it. Add a line to that list whenever a new interval is quoted outside a table.

## Layer 3 — superseded values in prose (`scripts/check_stale.py`)

Layer 2 reads tables and cannot see prose, so a quantity can be recomputed and its table
regenerated while a sentence or caption elsewhere still quotes the superseded value. This
layer holds every value the paper has superseded, with what replaced it and why it changed,
and fails the build if any of them survives anywhere in the source. A handful are quoted
deliberately, to show what a correction changed, so their permitted count is pinned too.

It cannot catch a number that was never right — only one that used to be right, which is
exactly the class Layer 2 misses.

## Continuous integration

`.github/workflows/reproduce.yml` performs the whole regeneration on a clean runner on
every push that touches the analysis code, after deleting the stored results so that the
checks cannot pass on them, and fails if any phase crashes even when the surviving numbers
match. Per-phase wall times are written to `results/runtimes.csv` on every run rather than
quoted in the paper, so the figure cannot go stale as phases are added.

The only components not reproducible without credentials are the two Materials Project
pulls, which require a free API key. Their cached outputs are committed and are labelled as
stored rather than regenerated wherever the CI run reuses them. See `results/PROVENANCE.md`.

## What none of this establishes

Reproducing a number is not the same as the number being right. All three layers test
consistency between code, results and manuscript, so they will pass on a statistical choice
that is internally consistent and still wrong. Every defect found in preparing this work
was found by reading the code against the manuscript, with all layers green. `REVIEWER.md`
separates the two questions, lists those defects as calibration, and names the definitions,
encodings, contrasts and the propensity model as the places most worth doubting.

## Live Materials Project re-pull (22 October 2026)

With `MP_API_KEY_FILE` set, `q2_known_baseline` and `q2b_completeness` query the live API
instead of using the committed pulls. A full run on 22 October 2026 (Windows, Python 3.12)
reproduced all 84 headline numbers with nothing stored. `check_paper.py` then reported 9
mismatches, all Materials Project drift since the 9–10 October pulls the paper reports:
`tab:mpcompl` boride/carbide/silicide n 2,365 → 2,367 (two r²SCAN entries added to MP),
and eight `tab:known` interval endpoints moved by ≤ 0.17 pp. No point estimate changed.
The committed pulls are the ones the paper reports; run `git checkout -- results` after a
live run so they are not replaced.

## Determinism across machines (checked 22 October 2026)

A full regeneration from a wiped `results/` on a second machine (Linux, Python 3.13,
pinned pandas 2.3.3 / numpy 2.2.6) reproduced all 84 headline quantities and all 655
printed values, and 51 of the 55 result CSVs byte-for-byte or to 1e-6. The four that
differ are `budget_model.csv` and `budget_repeated.csv` (early, superseded budget
evaluations; the paper reads `budget_baselines.csv` and `budget_strong*.csv`),
`fe_oxide_seed_stability.csv` (a 400-replicate seed probe; the paper reads
`percentile_convergence.csv`), and `selection_coverage.csv` (identical rows, different
sort order). No printed number reads any of them. `context_spread_ci.csv` reproduces only
because `run_all.py` runs `phase43_context_identify.py` twice with fixed seeds and
replicate counts (965 at 20261009, then 1,000 at 20261012); a single invocation gives
1,000 draws from one seed and moves the Table 3 endpoints by 1 to 7 meV/atom.

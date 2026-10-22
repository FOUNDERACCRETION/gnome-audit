# Reproducing every number in the paper

Two layers of checking run on every push; see `.github/workflows/reproduce.yml` for the
authoritative recipe. If you are reviewing rather than reproducing, read `REVIEWER.md`
instead — reproduction and correctness are different jobs.

## 1. Environment

```
python -m venv .venv
.venv/bin/pip install -r requirements.txt      # Windows: .venv\Scripts\pip
```

Python 3.10 or newer. pandas, numpy, matplotlib, and `requests` for the two Materials
Project pulls. **No scipy, no scikit-learn, no pymatgen.** The logistic regressions, the
cluster bootstraps and the sandwich estimator are written by hand in numpy, so no
dependency upgrade can silently change a number. `requirements.txt` pins the exact
versions every figure in the paper was produced with.

## 2. Data (not in the repo, ~175 MB)

Download both files from the GNoME public bucket into `data/`:

```
https://storage.googleapis.com/gdm_materials_discovery/gnome_data/stable_materials_r2scan.csv
https://storage.googleapis.com/gdm_materials_discovery/gnome_data/stable_materials_summary.csv
```

Verify the checksums, not just the row counts:

```
sha256sum data/*.csv
```

| File | Rows | Bytes | SHA-256 |
|---|---|---|---|
| `stable_materials_summary.csv` | 554,054 | 151,452,777 | `c37bffb9dd6c51d38d259a2d340ecf796194729f656d474828512390db1ff0c9` |
| `stable_materials_r2scan.csv` | 139,541 | 22,004,293 | `6a811395b065f43ae7911597405eec0ecfa88251e918f4c0ba8e4e99a75656fd` |

Both downloaded 9 October 2026. If a checksum differs, DeepMind has updated the file and
the numbers below will not match exactly; the paper states the counts and sums it used, so
a mismatch is informative rather than fatal. `MaterialId` is unique in both files, so no
duplicate resolution is required. Appendix E of the paper lists every column read, the
units, the missing-value coercion and the join rules.

## 3. Materials Project key (two phases only)

`q2_known_baseline.py` and `q2b_completeness.py` query the public API. Get a free key at
https://next-gen.materialsproject.org/api and export it:

```
export MP_API_KEY=...          # or: export MP_API_KEY_FILE=/path/to/file
```

Without a key those two phases are skipped and their cached outputs
(`results/known_baseline.csv`, `mp_completeness.csv`, `mp_gga_raw.csv`,
`mp_r2scan_raw.csv`) are reused. Everything that then reads them is **labelled
`[stored, not regenerated here]`** in the check output, so a run without a key cannot be
mistaken for a full one. The exact query — endpoint, the literal `thermo_types` values
`r2SCAN` and `GGA_GGA+U`, the `_fields` string and the paging — is in Appendix E.

## 4. Delete the stored results first

This is the step that makes the rest meaningful:

```
find results -type f \( -name '*.csv' -o -name '*.png' -o -name '*.json' \) \
     ! -name 'known_baseline.csv' ! -name 'mp_completeness.csv' \
     ! -name 'mp_gga_raw.csv' ! -name 'mp_r2scan_raw.csv' -delete
```

`results/` is committed, so without this the checks can pass on files you did not
generate. The CI workflow does exactly the above before every run.

## 5. Run

```
python run_all.py
```

This regenerates `results/` end to end, then checks two things:

- **84 headline quantities** against the values printed in the paper. Seven are
  claim-level rather than value-level — the budget ordering holding at every budget, the
  cross-fitted margin staying positive, the iron excess being smallest in oxides, the
  headline cell intervals staying disjoint, and so on. It deliberately does *not* assert
  that any interval excludes zero; that is a result, not an invariant.
- **655 printed table values**, via `scripts/check_paper.py`, which parses the numbers out
  of the LaTeX source and diffs them against the result files, interval endpoints
  included. A transcription error in a table fails the build.
- **Superseded values anywhere in the source**, via `scripts/check_stale.py`. The table
  checker cannot see prose, and the failure that has recurred most often in this paper is a
  recomputed quantity whose old value survives in a sentence or a caption. Every number that
  has changed is listed there with what replaced it; the check greps the whole source and
  fails on any hit, with an allowance for the few superseded values quoted deliberately as
  history.

A clean run ends with `ALL HEADLINE NUMBERS REPRODUCED` and `every parsed table value
matches the analysis output`. The driver **exits non-zero if any phase crashes even when
the surviving numbers match** — that failure mode has happened here, with a crashed phase
leaving stale outputs that the checks then passed.

Useful flags:

```
python run_all.py --check                       # only re-check existing results/
python run_all.py --skip phase33_percentile_conv,phase27_ame_ci    # drop the slow ones
```

## 6. Runtimes

`run_all.py` writes `results/runtimes.csv` with the wall time of every phase, the total,
the platform and the Python version. Read that rather than trusting a number in a README,
which goes stale the moment a phase is added. For orientation, the slowest phases are the
bootstrap and policy evaluations: `phase24_budget_strong` and `phase27_ame_ci` run two to
three minutes each, `phase16_excess_ci` about forty seconds at 100,000 replicates,
`phase33_percentile_conv` about a minute, `phase32` and `phase34` under a minute each.

Replicate counts are environment variables where it helps a reviewer: `NB` on
`phase16_excess_ci.py`, `phase32_bh_identity.py`, `phase34_converged_intervals.py` and
`phase35_coarse_cluster.py`, and `SEED` where a stability check wants a different draw.
Lowering `NB` makes a run fast and the intervals wrong; the committed defaults are the
ones the paper uses.

## 7. Rebuild the PDF

```
cd paper && pdflatex gnome_audit_paper.tex && pdflatex gnome_audit_paper.tex && pdflatex gnome_audit_paper.tex
```

Three passes, for the table of contents and the cross-references. A clean build reports no
errors, no overfull or underfull boxes and no undefined references. Figures resolve from
`paper/` first and then `../results/`, so a figure regenerated into `results/` will be
picked up only after you copy it into `paper/` or delete the stale copy there — a stale
figure has slipped through here before.

## 8. What the checks do and do not establish

`VALIDATION.md` records what each of the three layers asserts, why each of the 58 headline
checks is classified as `[impl]`, `[result]` or a reference comparison, and what none of
them can establish. Read it before concluding that a green run means the analysis is right.

## 9. If a number does not match

In order of likelihood:

1. You did not delete `results/` (§4), so you are checking files you did not make.
2. The input checksums differ (§2) — the upstream file changed.
3. A phase crashed. Read `results/<phase>.log`; the driver names the failing phase.
4. It is a real disagreement. That is the interesting case: please open an issue, and see
   `REVIEWER.md` for the parts of the analysis most worth doubting.

# Reproducing every number in the paper

Time on a 2-core laptop: about 40 minutes for the GNoME-only phases, plus 10–20 minutes
for the two Materials Project pulls.

## 1. Environment

```
python -m venv .venv
.venv/bin/pip install -r requirements.txt      # Windows: .venv\Scripts\pip
```

Python 3.10 or newer. No compiler, no scipy, no pymatgen.

## 2. Data (not in the repo, ~175 MB)

Download both files from the GNoME public bucket into `data/`:

```
https://storage.googleapis.com/gdm_materials_discovery/gnome_data/stable_materials_r2scan.csv
https://storage.googleapis.com/gdm_materials_discovery/gnome_data/stable_materials_summary.csv
```

Check:

```
python -c "import pandas as pd; r=pd.read_csv('data/stable_materials_r2scan.csv'); s=pd.read_csv('data/stable_materials_summary.csv', usecols=['MaterialId']); print(len(r), len(s))"
# expected: 139541 554054
```

If the row counts differ, DeepMind has updated the files and the numbers below will
not match exactly. The paper states the counts it used.

## 3. Materials Project key (only for §4.11 and §4.14)

Get a free key at https://next-gen.materialsproject.org/api and export it:

```
export MP_API_KEY=...          # or: set MP_API_KEY_FILE=/path/to/file
```

Without a key, `run_all.py` skips the two MP scripts and reports that §4.11 and the
completeness probe were not reproduced. Everything else runs.

## 4. Run

```
python run_all.py              # all phases, then the check
python run_all.py --check      # only compare existing results/ against the paper
python run_all.py --skip phase8_referee,q2b_completeness   # skip slow or network phases
```

`run_all.py` runs the scripts in dependency order, writes each one's stdout to
`results/<script>.log`, and finishes by comparing the headline quantities in
`results/` against the values printed in the paper (tolerance 0.1 percentage points
on rates, 0.02 on odds ratios; bootstrap phases are seeded so they should match
exactly, but floating-point order across BLAS builds can move the last digit).

## 5. What "reproduced" means here

The check compares these quantities, which are every number in the abstract:

| Quantity | Paper value | Source file |
|---|---|---|
| Global flip rate, MP reference | 13.04 % | `phase2_summary.json` |
| Global flip rate, MP+OQMD reference | 15.93 % | `phase2_summary.json` |
| Quaternary spread, top cell − bottom cell | 26.7 pp | `phase2_summary.json` |
| Top quaternary cell (oxide, 3d, no Ln) | 30.7 % | `stratified.csv` |
| Bottom quaternary cell (intermetallic, no 3d, no Ln) | 4.0 % | `stratified.csv` |
| Cluster-robust OR for has_U | 3.45 [3.13, 3.80] | `logit_cluster_robust.csv` |
| Budget: flips caught at 20 % recomputation, chemistry-ranked | 37.1 ± 0.9 % | `budget_repeated.csv` |
| Budget: random at 20 % | 20.1 % | `budget_repeated.csv` |
| Known-materials aggregate flip (MP, GGA-stable with r²SCAN) | 14.88 % | `known_baseline.csv` |

A re-run that matches all of these has reproduced the paper's central claims. The
secondary tables (per-element, band gap, Δ, hull-distance strata) are checked visually
against the paper's appendices; their scripts print them in the same layout.

## 6. If something does not match

Open an issue with the output of `python run_all.py --check` and
`pip freeze`. The most common cause is a changed upstream CSV (step 2).

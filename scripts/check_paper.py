#!/usr/bin/env python3
"""Verify that the numbers printed in the paper match the analysis output.

Reviewers of this kind of work have no way to tell whether a table was transcribed
correctly from the script that produced it. This closes that gap mechanically: it parses
the LaTeX source for the tables below, pairs each printed value with the corresponding
field in results/, and reports every disagreement. It runs in CI on every push.

    python scripts/check_paper.py [paper/gnome_audit_paper.tex]

Exit status 0 if every parsed value matches, 1 otherwise.
"""
import csv, re, sys, os
TEX = sys.argv[1] if len(sys.argv) > 1 else "paper/gnome_audit_paper.tex"
TOL = 0.051          # printed to one decimal, so half a unit in the last place
fails, checks = [], 0

def rows(f):
    with open(os.path.join("results", f)) as fh:
        return list(csv.DictReader(fh))

def body(label):
    """Return the row lines of the table carrying this \\label."""
    t = open(TEX).read()
    i = t.find("\\label{%s}" % label)
    if i < 0: return []
    seg = t[i:t.find("\\end{tab", i)]
    seg = seg[seg.find("\\midrule"):]
    return [l for l in seg.split("\\\\") if "&" in l]

def num(s):
    s = s.replace("{,}","").replace("\\%","").replace("\\phantom{0}","").replace("$","")
    s = s.replace("\\mathbf{","").replace("\\textbf{","").replace("}","")
    s = s.replace(",","").replace("+","").strip()
    m = re.search(r"-?\d+\.?\d*", s)
    return float(m.group()) if m else None

def cmp(label, got, want, tol=TOL):
    global checks
    checks += 1
    if got is None or want is None or abs(got-want) > tol:
        fails.append(f"{label}: paper {got}, results {want}")

SRC = open(TEX).read()

# ---- Table: within-anion excess (Section 4.3) -------------------------------
src = {r["anion"]: r for r in rows("within_anion_excess_ci.csv") if r["element"]=="Fe"}
ALL_EL = {(r["element"], r["anion"]): r for r in rows("within_anion_excess_ci.csv")}

# Within-anion intervals quoted in PROSE, not in a table. The table parser below cannot see
# these, and check_stale.py only catches values the paper has already superseded, so each
# prose-quoted interval is enumerated here by the sentence that carries it. Add a line here
# whenever a new element/anion interval is quoted outside tab:feexcess.
PROSE_INTERVALS = [
    ("Ni", "intermetallic", r"small excess in intermetallics \(\$\+([\d.]+)\$ \[\$\+([\d.]+)\$,\s*\$\+([\d.]+)\$\]\)"),
    ("Ni", "chalcogenide",  r"largest excesses are in chalcogenides \(\$\+([\d.]+)\$ \[\$\+([\d.]+)\$,\s*\$\+([\d.]+)\$\]\)"),
    ("Ni", "oxide",         r"chalcogenides \(\$[^)]*\]\) and oxides\s*\(\$\+([\d.]+)\$ \[\$\+([\d.]+)\$,\s*\$\+([\d.]+)\$\]\)"),
    ("Co", "oxide",         r"indistinguishable from zero\s*\(\$\+([\d.]+)\$ \[\$(-[\d.]+)\$,\s*\$\+([\d.]+)\$\]\)"),
    ("Fe", "oxide",         r"The oxide iron\s+excess, \$\+([\d.]+)\$ \[\$\+([\d.]+)\$,\s*\$\+([\d.]+)\$\]"),
]
for el, an, pat in PROSE_INTERVALS:
    m = re.search(pat, SRC)
    if m is None:
        checks += 1
        fails.append(f"prose interval for {el} in {an}: sentence not found (pattern changed?)")
    else:
        row = ALL_EL[(el, an)]
        cmp(f"prose {el} {an} excess",   float(m.group(1)), float(row["excess_pp"]))
        cmp(f"prose {el} {an} CI lower", float(m.group(2)), float(row["lo"]))
        cmp(f"prose {el} {an} CI upper", float(m.group(3)), float(row["hi"]))
for line in body("tab:feexcess"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 5 or "multicolumn" in line: continue
    name = c[0].replace("\\textbf{","").replace("}","").replace(" ($U$ applied)","").strip()
    name = {"boride/carbide/silicide":"boride/carbide/silicide"}.get(name,name)
    if name not in src: continue
    cmp(f"tab:feexcess {name} excess",   num(c[4]), float(src[name]["excess_pp"]))
    cmp(f"tab:feexcess {name} n",        num(c[1]), float(src[name]["n_with"]), 0.5)
    # the printed interval too, so a regenerated bootstrap cannot leave the table stale
    nums = re.findall(r"[-+]?\d+\.?\d*", c[4].replace("{,}","").replace("\\mathbf{","")
                                             .replace("\\textbf{","").replace("$",""))
    if len(nums) >= 3:
        cmp(f"tab:feexcess {name} CI lower", float(nums[1]), float(src[name]["lo"]))
        cmp(f"tab:feexcess {name} CI upper", float(nums[2]), float(src[name]["hi"]))

# ---- Table: cluster-robust joint test of equal interactions (Section 4.4) ---
jt = {r["element"]: r for r in rows("context_joint_test.csv")}
for line in body("tab:jointtest"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 4: continue
    el = c[0].replace("\\midrule", "").strip()
    if el not in jt: continue
    cmp(f"tab:jointtest {el} chi2", num(c[1]), float(jt[el]["wald_chi2"]), 0.011)
    cmp(f"tab:jointtest {el} BH threshold", num(c[3]), float(jt[el]["bh_threshold"]), 0.0011)

# ---- Table: context-dependence spreads and their intervals (Section 4.4) ----
cs = {r["element"]: r for r in rows("context_spread_ci.csv")}
cd = {}
for r_ in rows("context_dependence.csv"):
    cd.setdefault(r_["element"], {})[r_["anion"]] = float(r_["coef"])
CTX_COLS = ["oxide", "halide", "chalcogenide", "pnictide",
            "boride/carbide/silicide", "intermetallic"]
for line in body("tab:context"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 9: continue
    el = c[0].replace("\\midrule", "").strip()
    if el not in cs: continue
    # the six per-context effective coefficients, in the table's column order
    for i, an in enumerate(CTX_COLS):
        if an in cd.get(el, {}):
            cmp(f"tab:context {el} {an}", num(c[1 + i]), cd[el][an], 0.51)
    cmp(f"tab:context {el} spread", num(c[7]), float(cs[el]["spread_mev"]), 0.51)
    # strip \phantom{...} alignment padding first: its zeros are not data
    ci = re.sub(r"\\phantom\{[^}]*\}", "", c[8]).replace("$", "")
    nums = re.findall(r"[-+]?\d+\.?\d*", ci)
    if len(nums) >= 2:
        cmp(f"tab:context {el} spread CI lower", float(nums[0]), float(cs[el]["lo"]), 0.51)
        cmp(f"tab:context {el} spread CI upper", float(nums[1]), float(cs[el]["hi"]), 0.51)

# ---- Table: continuous hull-distance adjustment (Section 4.9) ---------------
hc = {(r["model"], r["feature"]): r for r in rows("hull_continuous.csv")}
HC_ROWS = [("$U$-list element present", "has_U"), ("other 3$d$ present", "has_3d_nonU"),
           ("lanthanide present", "has_Ln"), ("oxide", "is_oxide"),
           ("pnictide", "is_pnictide"), ("halide", "is_halide"),
           ("chalcogenide", "is_chalcogenide"),
           ("boride/carbide/silicide", "is_bcs"),
           ("McFadden pseudo-$R^2$", "mcfadden_r2")]
HC_MODELS = ((1, "no hull-distance term"),
             (2, "natural spline in PBE decomposition energy"))
for line in body("tab:hullcont"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 3: continue
    label = c[0].replace("\\midrule", "").strip()
    match = [k for lab, k in HC_ROWS if lab == label]
    if not match: continue
    key = match[0]
    for col, model in HC_MODELS:
        if col >= len(c): continue
        nums = re.findall(r"[-+]?\d+\.?\d*", c[col].replace("$", ""))
        if not nums: continue
        row = hc[(model, key)]
        cmp(f"tab:hullcont {label} [{model[:12]}] point",
            float(nums[0]), float(row["odds_ratio"]), 0.0051)
        if len(nums) >= 3 and row["or_lo"]:
            cmp(f"tab:hullcont {label} [{model[:12]}] lo", float(nums[1]), float(row["or_lo"]), 0.0051)
            cmp(f"tab:hullcont {label} [{model[:12]}] hi", float(nums[2]), float(row["or_hi"]), 0.0051)


# ---- Table: logistic regression (tab:logit): coef, OR, cluster-robust CI (variant B), SE inflation
lg = {r["feature"]: r for r in rows("logit_cluster_robust.csv")}
LG_ROWS = [("Intercept (baseline odds)", "intercept"), ("MP $U$-list element", "has_U"),
           ("Other 3$d$ (Sc, Ti, Cu, Zn)", "has_3d_nonU"), ("Lanthanide present", "has_Ln"),
           ("Oxide", "is_oxide"), ("Pnictide", "is_pnictide"), ("Boride/carbide/silicide", "is_bcs"),
           ("Chalcogenide", "is_chalc"), ("Halide", "is_halide"),
           ("Element count ($z$-scored)", "n_el_z"), ("$\\log$ sites ($z$-scored)", "lognsites_z")]
_seen = 0
for line in body("tab:logit"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 5: continue
    label = c[0].replace("\\midrule", "").strip()
    match = [k for lab, k in LG_ROWS if lab == label]
    if not match: continue
    r = lg[match[0]]; _seen += 1
    cmp(f"tab:logit {label} coef", num(c[1]), float(r["coef"]), 0.00051)
    cmp(f"tab:logit {label} OR", num(c[2]), float(r["OR"]), 0.00051)
    ci = re.findall(r"[-+]?\d+\.?\d*", c[3])
    cmp(f"tab:logit {label} lo", float(ci[0]), float(r["lo_B"]), 0.00051)
    cmp(f"tab:logit {label} hi", float(ci[1]), float(r["hi_B"]), 0.00051)
    cmp(f"tab:logit {label} SE inflation", num(c[4]), float(r["se_B"]) / float(r["se_naive"]), 0.0051)
if _seen != len(LG_ROWS):
    checks += 1; fails.append(f"tab:logit parsed {_seen} rows, expected {len(LG_ROWS)}")


# ---- Table: naive vs IPW flip rate by lanthanide-pooled quaternary cell (tab:ipw)
ipw = {(r["anion"], r["d3b"] == "True"): r for r in rows("ipw_cells.csv")}
_n = 0
for line in body("tab:ipw"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 4 or "spread" in c[0] or "ratio" in c[0]: continue
    lab = c[0].replace("\\midrule", "").strip()
    m = re.match(r"([a-z/]+), (3\$d\$|no 3\$d\$)", lab)
    if not m: continue
    r = ipw[(m.group(1), m.group(2) == "3$d$")]; _n += 1
    cmp(f"tab:ipw {lab} n", num(c[1]), float(r["n"]), 0.5)
    cmp(f"tab:ipw {lab} naive", num(c[2]), 100*float(r["naive"]))
    cmp(f"tab:ipw {lab} ipw", num(c[3]), 100*float(r["ipw"]))
if _n != 7: checks += 1; fails.append(f"tab:ipw parsed {_n} rows, expected 7")

# ---- Table: energy displacement by cell (tab:shift)
sh = {(r["anion"], r["d3"] == "True"): r for r in rows("shift_primary_cells.csv")}
SH_AN = {"oxide": "oxide", "bor./carb./sil.": "boride/carbide/silicide", "chalcogenide": "chalcogenide",
         "pnictide": "pnictide", "halide": "halide", "intermetallic": "intermetallic", "hydride": "hydride"}
_n = 0
for line in body("tab:shift"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 8: continue
    lab = c[0].replace("\\midrule", "").strip()
    d3 = "+ 3$d$" in lab
    an = SH_AN.get(lab.replace("+ 3$d$", "").strip())
    if an is None or (an, d3) not in sh: continue
    r = sh[(an, d3)]; _n += 1
    cmp(f"tab:shift {lab} n", num(c[1]), float(r["n"]), 0.5)
    for i, k in [(2, "mean"), (3, "median"), (4, "sd"), (5, "q1"), (6, "q3"), (7, "share_up")]:
        cmp(f"tab:shift {lab} {k}", num(c[i]), float(r[k]))
if _n < 8: checks += 1; fails.append(f"tab:shift parsed only {_n} rows")

# ---- Table: band-gap bins by 3d content (tab:gap)
bg = {r["gap_bin"]: r for r in rows("bandgap_bins_3d.csv")}
BG = {"0--0.5": "0-0.5", "0.5--1": "0.5-1", "1--3": "1-3", "$>3$": ">3"}
_n = 0
for line in body("tab:gap"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 6: continue
    lab = c[0].replace("\\midrule", "").strip()
    if lab not in BG: continue
    r = bg[BG[lab]]; _n += 1
    cmp(f"tab:gap {lab} n no3d", num(c[1]), float(r["n_no3d"]), 0.5)
    cmp(f"tab:gap {lab} flip no3d", num(c[2]), float(r["flip_no3d"]))
    cmp(f"tab:gap {lab} n 3d", num(c[3]), float(r["n_3d"]), 0.5)
    cmp(f"tab:gap {lab} flip 3d", num(c[4]), float(r["flip_3d"]))
    cmp(f"tab:gap {lab} excess", num(c[5]), float(r["diff_pp"]))
if _n != 4: checks += 1; fails.append(f"tab:gap parsed {_n} rows, expected 4")

# ---- Table: population flow (tab:popflow) against populations.csv
pop = {r["quantity"]: float(r["value"]) for r in rows("populations.csv")}
PF = [("Recomputation release", "r2scan_rows"), ("matched to the summary file", "matched"),
      ("matched and PBE-stable", "matched_pbe_stable"), ("Primary composition cells", "primary_cell_n")]
_n = 0
for line in body("tab:popflow"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 3: continue
    lab = c[0].replace("\\quad", "").replace("\\midrule", "").strip()
    for lab2, k in PF:
        if lab == lab2:
            _n += 1; cmp(f"tab:popflow {lab}", num(c[1]), pop[k], 0.5)
if _n != len(PF): checks += 1; fails.append(f"tab:popflow matched {_n} rows, expected {len(PF)}")


# ---- Table: flip rate against threshold (tab:appthreshold)
th = {float(r["threshold_eV"]): r for r in rows("threshold_curves.csv")}
TH_COLS = ["overall", "boride/carbide/silicide", "chalcogenide", "halide", "hydride", "intermetallic", "oxide"]
_n = 0
for line in body("tab:appthreshold"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 8: continue
    try: tau = float(c[0].replace("\\midrule", "").strip())
    except ValueError: continue
    if tau not in th: continue
    r = th[tau]; _n += 1
    for i, k in enumerate(TH_COLS):
        cmp(f"tab:appthreshold tau={tau} {k}", num(c[i+1]), 100*float(r[k]))
if _n != 7: checks += 1; fails.append(f"tab:appthreshold parsed {_n} rows, expected 7")

# ---- Table: headline cells under both reference sets (Section 4.16) ---------
ar = {(r["reference"], r["cell"]): r for r in rows("allref_cells.csv")}
AR_ROWS = [("oxide, no Ln, 3$d$ (4-el)", "oxide, no Ln, 3d"),
           ("intermetallic, no Ln/3$d$ (4-el)", "intermetallic, no Ln/3d"),
           ("spread (pp)", "spread (oxide - intermetallic, pp)"),
           ("ratio", "ratio (oxide / intermetallic)")]
for line in body("tab:allcells"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 4: continue
    label = c[0].replace("\\phantom{0}", "").replace("\\midrule", "").strip()
    match = [k for lab, k in AR_ROWS if lab == label]
    if not match: continue
    key = match[0]
    for col, ref in ((2, "MP"), (3, "All")):
        txt = c[col].replace("{,}", "").replace("\\%", "").replace("\\phantom{0}", "")
        nums = re.findall(r"[-+]?\d+\.?\d*", txt.replace("$", ""))
        if not nums: continue
        row = ar[(ref, key)]
        cmp(f"tab:allcells {label} {ref} point", float(nums[0]), float(row["rate_pct"]))
        # the MP-reference cell rates are printed as point estimates only, by design
        if len(nums) >= 3:
            cmp(f"tab:allcells {label} {ref} CI lower", float(nums[1]), float(row["lo"]))
            cmp(f"tab:allcells {label} {ref} CI upper", float(nums[2]), float(row["hi"]))

# ---- Table: matched cells against known materials (Section 4.13) ------------
kc = list(rows("known_cells_clustered_release.csv"))
def _cellkey(name):
    """'4-el oxide, 3$d$' -> (n_el, anion, has_Ln, has_3d) as the CSV spells them."""
    n = name.split("-el")[0].strip()
    parts = [x.strip() for x in name.split("-el", 1)[1].split(",")]
    rest, flags = parts[0], parts[1:]
    # the row label spells the flags as ", Ln", ", 3$d$", ", Ln, 3$d$" or ", no Ln/3$d$";
    # the last of those contains the substring "3$d$" and must not be read as the 3d flag
    if any(f.startswith("no Ln") for f in flags):
        ln = d3 = False
    else:
        ln = any(f == "Ln" for f in flags)
        d3 = any(f == "3$d$" for f in flags)
    for a, short in (("boride/carbide/silicide", "bor./carb./sil."), ("chalcogenide", "chalcogenide"),
                     ("intermetallic", "intermetallic"), ("oxide", "oxide"), ("halide", "halide"),
                     ("pnictide", "pnictide"), ("hydride", "hydride")):
        if short in rest:
            return (n, a, str(ln), str(d3))
    return None
for line in body("tab:known"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 7 or "multicolumn" in line: continue
    k = _cellkey(c[0])
    if k is None: continue
    src = [r for r in kc if (r["n_elements"], r["anion_class"], r["has_Ln"], r["has_3d"]) == k]
    if not src: continue
    r0 = src[0]
    cmp(f"tab:known {c[0]} n_G",   num(c[1]), float(r0["n_G"]), 0.5)
    cmp(f"tab:known {c[0]} delta", num(c[5]), float(r0["delta"]))
    nums = re.findall(r"[-+]?\d+\.?\d*", c[6].replace("{,}", "").replace("$", ""))
    if len(nums) >= 2:
        cmp(f"tab:known {c[0]} CI lower", float(nums[0]), float(r0["lo"]))
        cmp(f"tab:known {c[0]} CI upper", float(nums[1]), float(r0["hi"]))

# ---- Table: element x anion context (Section 4.4) ---------------------------
ctx = {}
for r in rows("context_dependence.csv"):
    ctx.setdefault(r["element"], {})[r["anion"]] = float(r["coef"])
order = ["oxide","halide","chalcogenide","pnictide","boride/carbide/silicide","intermetallic"]
for line in body("tab:context"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 8: continue
    el = c[0].strip()
    if el not in ctx: continue
    for k, anion in enumerate(order, start=1):
        v = num(c[k])
        if v is None or anion not in ctx[el]: continue
        cmp(f"tab:context {el}/{anion}", v, ctx[el][anion], 0.51)

# ---- Table: budget policies (Section 4.14) ----------------------------------
bud = {round(float(r["budget"]),2): r for r in rows("budget_baselines.csv")}
for line in body("tab:budget"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 5: continue
    b = num(c[0])
    if b is None: continue
    key = round(b/100, 2)
    if key not in bud: continue
    for k, col in zip(range(1,5), ["random","chemistry","hull","combined"]):
        cmp(f"tab:budget {b:.0f}%/{col}", num(c[k]), float(bud[key][col])*100)

# ---- Table: MP hull completeness (Section 4.16) -----------------------------
comp = {}
for r in rows("mp_completeness.csv"):
    comp.setdefault(r["anion"], []).append(r)
for line in body("tab:mpcompl"):
    c = [x.strip() for x in line.split("&")]
    if len(c) < 4: continue
    a = c[0].strip()
    if a not in comp: continue
    v = comp[a]
    cmp(f"tab:mpcompl {a} n", num(c[1]), float(len(v)), 0.5)
    flips = sum(int(x["flip"]) for x in v)
    cmp(f"tab:mpcompl {a} flip", num(c[3]), flips/len(v)*100)

# ---- Appendix: per-element -------------------------------------------------
pe = {r["element"]: r for r in rows("per_element_clustered.csv")}
seen = 0
t = open(TEX).read()
i = t.find("\\section{Per-element results}")
if i > 0:
    seg = t[i:t.find("\\end{longtable}", i)]
    for m in re.finditer(r"([A-Z][a-z]?)(?:\$\^\{\*\}\$)?\s*&\s*((?:\d|\{,\})+)\s*&\s*([\d.]+)", seg):
        el, n, f = m.group(1), int(m.group(2).replace("{,}", "")), float(m.group(3))
        if el not in pe: continue
        seen += 1
        cmp(f"app:element {el} n", float(n), float(pe[el]["n"]), 0.5)
        cmp(f"app:element {el} flip", f, float(pe[el]["flip"])*100)
    if seen != len(pe):
        fails.append(f"app:element lists {seen} elements, results has {len(pe)}")
        checks += 1

# Widening statistics quoted in prose (Section 3): the median and range of the cluster/naive
# interval-width ratio over the primary composition cells.
_w = sorted(float(r["widening"]) for r in rows("cell_ci_converged.csv"))
_m = re.search(r"observed widening relative to naive intervals has median ([\d.]+) and ranges from\s+([\d.]+) to\s+([\d.]+)", SRC)
if _m is None:
    checks += 1; fails.append("prose widening sentence not found (pattern changed?)")
else:
    _med = _w[len(_w)//2] if len(_w) % 2 else 0.5*(_w[len(_w)//2-1]+_w[len(_w)//2])
    cmp("prose widening median", float(_m.group(1)), _med, 0.0051)
    cmp("prose widening min",    float(_m.group(2)), _w[0], 0.0051)
    cmp("prose widening max",    float(_m.group(3)), _w[-1], 0.0051)
    cmp("prose widening cell count", 40, len(_w), 0.5)

print(f"checked {checks} printed values against results/")
if fails:
    print(f"\n{len(fails)} MISMATCH(ES):")
    for f in fails: print("  " + f)
    sys.exit(1)
print("every parsed table value matches the analysis output")

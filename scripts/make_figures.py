#!/usr/bin/env python3
"""Regenerate the four manuscript figures that are drawn from results/ tables.

    python scripts/make_figures.py

fig1_flip_by_chemistry.png  from results/cell_ci_converged.csv and results/phase2_summary.json
fig2_budget.png             from results/budget_baselines.csv
fig3_gnome_vs_known.png     from results/known_cells_clustered_release.csv
fig4_per_element.png        from results/per_element_clustered.csv and phase2_summary.json

All four are read-only views of tables the paper already prints, so no number appears
here that is not also checked by check_paper.py. Figure sizes are the physical print
sizes used by the manuscript (text width 6.5 in), so fonts are not scaled down on insertion.
"""
import csv, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7})

OUT = "results"

def rows(f):
    with open(os.path.join(OUT, f)) as fh:
        return list(csv.DictReader(fh))

def tf(s):
    return str(s).strip().lower() in ("true", "1", "yes")

# ---------------------------------------------------------------- fig1
cells = [r for r in rows("cell_ci_converged.csv") if int(r["n_el"]) == 4]
cells.sort(key=lambda r: float(r["rate"]))
glob = 100 * json.load(open(os.path.join(OUT, "phase2_summary.json")))["global_flip_mp"]

def label(r):
    s = r["anion"].replace("boride/carbide/silicide", "B/C/Si")
    if tf(r["Ln"]): s += " +Ln"
    if tf(r["d3"]): s += " +3d"
    return s

fig, ax = plt.subplots(figsize=(3.4, 6.6), dpi=300)
y = range(len(cells))
rate = [100 * float(r["rate"]) for r in cells]
ax.barh(list(y), rate, color="#4C72B0", height=0.8)
for i, r in enumerate(cells):
    nlo, nhi = 100 * float(r["naive_lo"]), 100 * float(r["naive_hi"])
    clo, chi = 100 * float(r["clust_lo"]), 100 * float(r["clust_hi"])
    ax.plot([nlo, nhi], [i, i], color="0.55", lw=5, solid_capstyle="butt", zorder=3)
    ax.errorbar(100 * float(r["rate"]), i, xerr=[[100 * float(r["rate"]) - clo], [chi - 100 * float(r["rate"])]],
                fmt="none", ecolor="black", elinewidth=1.3, capsize=3, zorder=4)
ax.axvline(glob, color="#C44E52", ls="--", lw=1.2, zorder=2,
           label=f"all element counts, {glob:.2f}%")
ax.legend(loc="lower right", frameon=True)
ax.set_yticks(list(y))
ax.set_yticklabels([label(r) for r in cells])
ax.set_xlabel("PBE→r²SCAN flip rate (%)")
ax.set_xlim(0, 36)
ax.set_ylim(-0.8, len(cells) - 0.2)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_flip_by_chemistry.png"))
plt.close(fig)

# ---------------------------------------------------------------- fig3
known = rows("known_cells_clustered_release.csv")

def short(r):
    s = f"{r['n_elements']}-el {r['anion_class']}"
    s = s.replace("boride/carbide/silicide", "bor./carb./sil.")
    if tf(r["has_Ln"]): s += ", Ln"
    if tf(r["has_3d"]): s += ", 3d"
    return s

fig, ax = plt.subplots(figsize=(4.6, 4.4), dpi=300)
ax.plot([0, 35], [0, 35], ls="--", color="0.5", lw=1)
for r in known:
    xk, yg = 100 * float(r["flip_known"]), 100 * float(r["flip_G"])
    ax.scatter(xk, yg, s=28, color="#4C72B0", edgecolor="white", lw=0.5, zorder=3)
    if tf(r["excl0"]):
        # the crowded cells at the lower left get explicit offsets so labels do not collide
        s = short(r)
        if s == "4-el oxide, 3d":
            ax.annotate(s, (xk, yg), xytext=(-5, 0), textcoords="offset points", fontsize=7.5, ha="right", va="center")
        elif s == "3-el chalcogenide, 3d":
            ax.annotate(s, (xk, yg), xytext=(4, 9), textcoords="offset points", fontsize=7.5,
                        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6))
        elif s == "3-el bor./carb./sil., 3d":
            ax.annotate(s, (xk, yg), xytext=(6, -11), textcoords="offset points", fontsize=7.5,
                        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6))
        elif s == "3-el pnictide":
            ax.annotate(s, (xk, yg), xytext=(-18, 10), textcoords="offset points", fontsize=7.5, ha="right",
                        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6))
        elif s == "3-el chalcogenide":
            ax.annotate(s, (xk, yg), xytext=(-18, -14), textcoords="offset points", fontsize=7.5, ha="right",
                        arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6))
        elif xk < 8.5:
            ax.annotate(s, (xk, yg), xytext=(-4, -9), textcoords="offset points", fontsize=7.5, ha="right")
        else:
            ax.annotate(s, (xk, yg), xytext=(5, 3), textcoords="offset points", fontsize=7.5)
ax.set_xlabel("previously catalogued materials (Materials Project), flip rate (%)")
ax.set_ylabel("GNoME flip rate (%)")
ax.set_xlim(0, 35)
ax.set_ylim(0, 36)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3_gnome_vs_known.png"))
plt.close(fig)
print("wrote fig1_flip_by_chemistry.png and fig3_gnome_vs_known.png")

# ---------------------------------------------------------------- fig2
bb = rows("budget_baselines.csv")
b = [0.0] + [100 * float(r["budget"]) for r in bb]
fig, ax = plt.subplots(figsize=(4.6, 3.5), dpi=300)
pol = [("combined", "both", "#4C72B0", "-", "o"), ("hull", "hull distance", "#55A868", "-", "s"),
       ("chemistry", "composition", "#C44E52", "-", "^"), ("random", "uniform random", "0.45", "--", "D")]
for key, lab, col, ls, mk in pol:
    m = [0.0] + [100 * float(r[key]) for r in bb]
    s = [0.0] + [100 * float(r[key + "_sd"]) for r in bb]
    ax.plot(b, m, ls=ls, marker=mk, color=col, lw=1.5, ms=4, label=lab, zorder=3)
    ax.fill_between(b, [x - y for x, y in zip(m, s)], [x + y for x, y in zip(m, s)], color=col, alpha=0.18, lw=0)
ax.set_xlabel("r²SCAN recomputation budget (% of materials)")
ax.set_ylabel("disagreements surfaced (%)")
ax.set_xlim(0, 52); ax.set_ylim(0, 90)
ax.grid(alpha=0.3); ax.legend(loc="upper left", frameon=True)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_budget.png")); plt.close(fig)

# ---------------------------------------------------------------- fig4
pe = [r for r in rows("per_element_clustered.csv") if int(float(r["n"])) >= 300]
pe.sort(key=lambda r: -float(r["flip"]))
top = pe[:30]
ULIST = {"V", "Cr", "Mn", "Fe", "Co", "Ni", "W", "Mo"}
fig, ax = plt.subplots(figsize=(3.0, 6.6), dpi=300)
ys = list(range(len(top)))[::-1]
for y, r in zip(ys, top):
    f, lo, hi = 100 * float(r["flip"]), 100 * float(r["clust_lo"]), 100 * float(r["clust_hi"])
    ax.barh(y, f, color="#C44E52" if r["element"] in ULIST else "#4C72B0", height=0.8)
    ax.errorbar(f, y, xerr=[[f - lo], [hi - f]], fmt="none", ecolor="black", elinewidth=1.2, capsize=3, zorder=4)
ax.axvline(glob, color="0.5", ls="--", lw=1.2, zorder=2)
ax.set_yticks(ys); ax.set_yticklabels([f"{r['element']} (n={int(float(r['n'])):,})" for r in top])
ax.set_xlabel("flip rate (%)")
ax.set_xlim(0, 48); ax.set_ylim(-0.8, len(top) - 0.2)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_per_element.png")); plt.close(fig)
print("wrote fig2_budget.png and fig4_per_element.png")

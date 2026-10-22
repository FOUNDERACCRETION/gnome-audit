"""Fig 2 redrawn: four ranking policies, from results/budget_baselines.csv (means over
20 cluster-wise splits), replacing the single-split chemistry-vs-random curve."""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open("results/budget_baselines.csv")))
b=[float(r["budget"])*100 for r in rows]
ser=[("Both",        [float(r["combined"])*100 for r in rows], "#1b4965", "o", "-"),
     ("Hull distance",[float(r["hull"])*100     for r in rows], "#5fa8d3", "s", "-"),
     ("Composition", [float(r["chemistry"])*100 for r in rows], "#c1666b", "^", "-"),
     ("Uniform random",[float(r["random"])*100  for r in rows], "#8d8d8d", None, "--")]
fig,ax=plt.subplots(figsize=(5.6,4.0))
for lab,v,c,mk,ls in ser:
    ax.plot(b,v,ls,color=c,marker=mk,ms=5,lw=1.8,label=lab)
    if lab in ("Composition","Hull distance","Both"):
        sd=[float(r[{"Composition":"chemistry_sd","Hull distance":"hull_sd","Both":"combined_sd"}[lab]])*100 for r in rows]
        ax.fill_between(b,[x-s for x,s in zip(v,sd)],[x+s for x,s in zip(v,sd)],color=c,alpha=.15,lw=0)
ax.set_xlabel("r$^2$SCAN recomputation budget (% of candidates)")
ax.set_ylabel("Disagreements surfaced (%)")
ax.set_title("Hull distance beats composition; together they beat both",fontsize=9.5)
ax.legend(fontsize=8,frameon=False); ax.grid(alpha=.3); ax.set_xlim(0,52); ax.set_ylim(0,90)
fig.tight_layout(); fig.savefig("results/fig2_budget.png",dpi=200)
print("-> results/fig2_budget.png")

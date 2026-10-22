"""Regenerate the band-gap and logistic-regression figures so they match the tables.
The previous versions predated the Section 4.8 rebinning and the relabelling of the
U-list regressor, and plotted naive intervals the text tells the reader not to use.
"""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
R="results"

# ---- Fig 3: flip rate by band-gap bin, split on 3d, matching tab:gap ----------
rows=list(csv.DictReader(open(f"{R}/bandgap_bins_3d.csv")))
labs=[r["gap_bin"] for r in rows]
a=[float(r["flip_no3d"]) for r in rows]; b=[float(r["flip_3d"]) for r in rows]
na=[int(r["n_no3d"]) for r in rows];     nb=[int(r["n_3d"]) for r in rows]
x=np.arange(len(labs)); w=0.38
fig,ax=plt.subplots(figsize=(5.8,3.9))
r1=ax.bar(x-w/2,a,w,label="no 3$d$ element",color="#8aa8c0")
r2=ax.bar(x+w/2,b,w,label="3$d$ element present",color="#b5484f")
for xi,v,n in zip(x-w/2,a,na): ax.text(xi,v+0.5,f"{v:.1f}%\n$n$={n:,}",ha="center",fontsize=6.3)
for xi,v,n in zip(x+w/2,b,nb): ax.text(xi,v+0.5,f"{v:.1f}%\n$n$={n:,}",ha="center",fontsize=6.3)
ax.set_xticks(x); ax.set_xticklabels([f"{l} eV" for l in labs])
ax.set_xlabel("PBE band gap"); ax.set_ylabel("Flip rate (%)")
ax.set_title("Flip rate by PBE band-gap bin and 3$d$ content",fontsize=9)
ax.set_ylim(0,38); ax.legend(fontsize=8,frameon=False); ax.grid(axis="y",alpha=.3)
fig.tight_layout(); fig.savefig(f"{R}/fig5_bandgap.png",dpi=200); plt.close(fig)
print("-> fig5_bandgap.png")

# ---- Fig 4: odds ratios with CLUSTER-ROBUST intervals and the table's labels ---
NAME={"has_U":"MP $U$-list element","has_3d_nonU":"Other 3$d$ (Sc, Ti, Cu, Zn)",
      "has_Ln":"Lanthanide present","is_oxide":"Oxide","is_pnictide":"Pnictide",
      "is_bcs":"Boride/carbide/silicide","is_chalc":"Chalcogenide","is_halide":"Halide",
      "n_el_z":"Element count ($z$)","lognsites_z":"$\\log$ sites ($z$)"}
rows=[r for r in csv.DictReader(open(f"{R}/logit_cluster_robust.csv"))
      if r["feature"] in NAME]
rows.sort(key=lambda r: float(r["OR"]))
y=np.arange(len(rows))
orr=[float(r["OR"]) for r in rows]
lo=[float(r["lo_B"]) for r in rows]; hi=[float(r["hi_B"]) for r in rows]
fig,ax=plt.subplots(figsize=(5.8,4.0))
ax.errorbar(orr,y,xerr=[np.array(orr)-np.array(lo),np.array(hi)-np.array(orr)],
            fmt="o",ms=4.5,lw=1.3,capsize=3,color="#1b4965")
ax.axvline(1,color="#999",ls="--",lw=1)
ax.set_yticks(y); ax.set_yticklabels([NAME[r["feature"]] for r in rows],fontsize=8)
ax.set_xscale("log"); ax.set_xlabel("Odds ratio (log scale), cluster-robust 95 % CI")
ax.set_title("Cluster-robust odds ratios; reference is a plain intermetallic",fontsize=9)
ax.grid(axis="x",alpha=.3)
fig.tight_layout(); fig.savefig(f"{R}/fig6_logit.png",dpi=200); plt.close(fig)
print("-> fig6_logit.png")

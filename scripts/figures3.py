import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
R="results"
# Fig 1 redone with CLUSTER-ROBUST intervals
rows=[r for r in csv.DictReader(open(f"{R}/cluster_bootstrap.csv")) if r["n_el"]=="4"]
rows.sort(key=lambda r: float(r["rate"]))
lab=[f"{r['anion']}{' +Ln' if r['Ln']=='True' else ''}{' +3d' if r['d3']=='True' else ''}"
     for r in rows]
v=np.array([float(r["rate"]) for r in rows])*100
lo=np.array([float(r["clust_lo"]) for r in rows])*100
hi=np.array([float(r["clust_hi"]) for r in rows])*100
nlo=np.array([float(r["naive_lo"]) for r in rows])*100
nhi=np.array([float(r["naive_hi"]) for r in rows])*100
fig,ax=plt.subplots(figsize=(8,9)); y=np.arange(len(rows))
ax.barh(y,v,color="#4C72B0",zorder=2)
ax.errorbar(v,y,xerr=[v-nlo,nhi-v],fmt="none",ecolor="#B0B0B0",elinewidth=3,capsize=0,zorder=3)
ax.errorbar(v,y,xerr=[v-lo,hi-v],fmt="none",ecolor="black",elinewidth=1.1,capsize=2.5,zorder=4)
ax.axvline(13.04,color="#C44E52",ls="--",lw=1.2,zorder=1)
ax.set_yticks(y); ax.set_yticklabels(lab,fontsize=8)
ax.set_xlabel("PBE→r²SCAN flip rate (%)")
ax.set_title("Quaternary GNoME materials: flip rate by chemistry\n"
  "black = cluster-robust 95% CI, grey = naive binomial (too narrow)",fontsize=10)
fig.tight_layout(); fig.savefig(f"{R}/fig1_flip_by_chemistry.png",dpi=200)

# Fig 4 redone with cluster-robust intervals
pe=[r for r in csv.DictReader(open(f"{R}/per_element_clustered.csv"))]
pe.sort(key=lambda r:-float(r["flip"])); pe=pe[:30][::-1]
v=np.array([float(r["flip"]) for r in pe])*100
lo=np.array([float(r["clust_lo"]) for r in pe])*100
hi=np.array([float(r["clust_hi"]) for r in pe])*100
U=set("Co Cr Fe Mn Mo Ni V W".split())
fig,ax=plt.subplots(figsize=(7,9)); y=np.arange(len(pe))
ax.barh(y,v,xerr=[v-lo,hi-v],
        color=["#C44E52" if r["element"] in U else "#4C72B0" for r in pe],
        error_kw=dict(lw=1,capsize=2))
ax.axvline(13.04,color="grey",ls="--",lw=1.2)
ax.set_yticks(y); ax.set_yticklabels([f"{r['element']} (n={r['n']})" for r in pe],fontsize=8)
ax.set_xlabel("flip rate (%)")
ax.set_title("Per-element flip rate, top 30 (n≥300)\n"
  "cluster-robust 95% CI; red = elements receiving MP2020 $U$",fontsize=10)
fig.tight_layout(); fig.savefig(f"{R}/fig4_per_element.png",dpi=200)
print("regenerated fig1, fig4 with cluster-robust intervals")

import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
R="results"
def rd(f): return list(csv.DictReader(open(f"{R}/{f}")))

# Fig 1: flip rate by cell, n_elements=4, with CIs
s=[r for r in rd("stratified.csv") if r["n_elements"]=="4"]
s.sort(key=lambda r: float(r["flip_mp"]))
lab=[f"{r['anion_class']}{' +Ln' if r['has_Ln']=='True' else ''}"
     f"{' +3d' if r['has_3d']=='True' else ''}" for r in s]
v=np.array([float(r["flip_mp"]) for r in s])*100
lo=np.array([float(r["lo"]) for r in s])*100; hi=np.array([float(r["hi"]) for r in s])*100
fig,ax=plt.subplots(figsize=(8,9))
y=np.arange(len(s))
ax.barh(y,v,xerr=[v-lo,hi-v],color="#4C72B0",error_kw=dict(lw=1,capsize=2))
ax.axvline(13.04,color="#C44E52",ls="--",lw=1.2,label="global 13.0%")
ax.set_yticks(y); ax.set_yticklabels(lab,fontsize=8)
ax.set_xlabel("PBE→r²SCAN flip rate (%)")
ax.set_title("Quaternary GNoME materials: stability-verdict flip rate by chemistry\n"
             "(disjoint cells, n≥500, bootstrap 95% CI)",fontsize=10)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{R}/fig1_flip_by_chemistry.png",dpi=200)

# Fig 2: budget curve
b=rd("budget_model.csv")
x=[float(r["budget"])*100 for r in b]
c=[float(r["chem"])*100 for r in b]; rn=[float(r["random"])*100 for r in b]
fig,ax=plt.subplots(figsize=(6,4.5))
ax.plot([0]+x,[0]+c,"o-",color="#4C72B0",label="chemistry-ranked")
ax.plot([0]+x,[0]+rn,"s--",color="#937860",label="uniform random")
ax.set_xlabel("r²SCAN recomputation budget (% of materials)")
ax.set_ylabel("wrong stability labels caught (%)")
ax.set_title("Verification budget: chemistry-aware vs uniform allocation",fontsize=10)
ax.legend(fontsize=9); ax.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{R}/fig2_budget.png",dpi=200)

# Fig 3: GNoME vs known, matched cells
k={(r["n_elements"],r["anion_class"],r["has_Ln"],r["has_3d"]):r for r in rd("known_baseline.csv")}
pairs=[]
for r in rd("stratified.csv"):
    key=(r["n_elements"],r["anion_class"],r["has_Ln"],r["has_3d"])
    if key in k:
        pairs.append((float(r["flip_mp"])*100,float(k[key]["flip_known"])*100,
                      f"{r['n_elements']}el {r['anion_class'][:12]}"
                      f"{'+Ln' if r['has_Ln']=='True' else ''}{'+3d' if r['has_3d']=='True' else ''}"))
fig,ax=plt.subplots(figsize=(6.5,6))
g=[p[0] for p in pairs]; kk=[p[1] for p in pairs]
ax.scatter(kk,g,c=["#C44E52" if a>b else "#55A868" for a,b in zip(g,kk)],s=45,zorder=3)
m=max(max(g),max(kk))*1.1
ax.plot([0,m],[0,m],color="grey",lw=1,ls="--",zorder=1)
for a,b,l in pairs:
    if abs(a-b)>5: ax.annotate(l,(b,a),fontsize=6.5,xytext=(3,3),textcoords="offset points")
ax.set_xlabel("known materials (Materials Project) flip rate (%)")
ax.set_ylabel("GNoME flip rate (%)")
ax.set_title("Matched chemistry cells: AI-discovered vs known\n"
             "red = GNoME worse, green = GNoME better",fontsize=10)
ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(f"{R}/fig3_gnome_vs_known.png",dpi=200)
print("wrote 3 figures")

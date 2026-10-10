import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
R="results"
def rd(f): return list(csv.DictReader(open(f"{R}/{f}")))

# Fig 4: per-element, top 30
pe=rd("per_element.csv")
pe=sorted(pe,key=lambda r:-float(r["flip"]))[:30][::-1]
v=np.array([float(r["flip"]) for r in pe])*100
lo=np.array([float(r["lo"]) for r in pe])*100; hi=np.array([float(r["hi"]) for r in pe])*100
fig,ax=plt.subplots(figsize=(7,9)); y=np.arange(len(pe))
cols=["#C44E52" if r["element"] in ("Fe","Co","Cr","Mn","Mo","Ni","V","W") else "#4C72B0" for r in pe]
ax.barh(y,v,xerr=[v-lo,hi-v],color=cols,error_kw=dict(lw=1,capsize=2))
ax.axvline(13.04,color="grey",ls="--",lw=1.2)
ax.set_yticks(y); ax.set_yticklabels([f"{r['element']} (n={r['n']})" for r in pe],fontsize=8)
ax.set_xlabel("flip rate (%)")
ax.set_title("Per-element flip rate, top 30 (n$\\geq$300)\nred = elements receiving MP2020 $U$ in oxides/fluorides",fontsize=10)
fig.tight_layout(); fig.savefig(f"{R}/fig4_per_element.png",dpi=200)

# Fig 5: band gap bins
gb=rd("bandgap_bins.csv")
lab=[r["gapbin"] for r in gb]; v=[float(r["flip"])*100 for r in gb]; n=[int(r["n"]) for r in gb]
fig,ax=plt.subplots(figsize=(6.5,4.2))
ax.bar(range(len(gb)),v,color="#4C72B0")
for i,(a,b) in enumerate(zip(v,n)): ax.text(i,a+.3,f"n={b}",ha="center",fontsize=7)
ax.axhline(13.04,color="#C44E52",ls="--",lw=1.2,label="global 13.0%")
ax.set_xticks(range(len(gb))); ax.set_xticklabels(["metal\n(0)","0–0.5","0.5–1.5","1.5–3","> 3"],fontsize=9)
ax.set_xlabel("PBE band gap (eV)"); ax.set_ylabel("flip rate (%)")
ax.set_title("Flip rate against band gap: non-monotonic",fontsize=10)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{R}/fig5_bandgap.png",dpi=200)

# Fig 6: logistic odds ratios
lg=[r for r in rd("logistic_coefficients.csv") if r["feature"]!="intercept"][::-1]
o=np.array([float(r["odds_ratio"]) for r in lg])
l=np.array([float(r["or_lo"]) for r in lg]); h=np.array([float(r["or_hi"]) for r in lg])
fig,ax=plt.subplots(figsize=(6.5,4.6)); y=np.arange(len(lg))
ax.errorbar(o,y,xerr=[o-l,h-o],fmt="o",color="#4C72B0",capsize=3)
ax.axvline(1,color="grey",ls="--",lw=1)
ax.set_yticks(y); ax.set_yticklabels([r["feature"] for r in lg],fontsize=8)
ax.set_xlabel("odds ratio for a flipped verdict (95% CI)")
ax.set_title("Logistic model of flip probability from composition",fontsize=10)
ax.set_xscale("log"); fig.tight_layout(); fig.savefig(f"{R}/fig6_logit.png",dpi=200)
print("wrote fig4, fig5, fig6")

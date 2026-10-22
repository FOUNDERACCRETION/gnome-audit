"""Items 37/38: direct standardisation of the GNoME/MP aggregate comparison over the
matched cells actually used in Section 4.13 (n_G>=500, n_K>=100), plus a LaTeX
supplementary table of all of them. Also writes the AME intervals of item 15."""
import numpy as np, pandas as pd
j=pd.read_csv("results/known_all_cells.csv")
full=j.copy()
j=j[(j.n_G>=500)&(j.n>=100)].copy().sort_values("delta",ascending=False)
print("cells in the Section 4.13 set: %d  (worse %d, better %d, intervals excluding zero %d)"
      %(len(j),int((j.delta>0).sum()),int((j.delta<0).sum()),int(j.excl0.sum())))
wG=j.n_G/j.n_G.sum(); wK=j.n/j.n.sum()
cg=(j.flip_G*wG).sum()*100; ck=(j.flip_known*wK).sum()*100
sg=(j.flip_G*wK).sum()*100; sk=(j.flip_known*wG).sum()*100
print("crude GNoME over these cells        %.2f%%"%cg)
print("crude Known over these cells        %.2f%%"%ck)
print("GNoME standardised to Known mix     %.2f%%   (difference %+.2f pp)"%(sg,sg-ck))
print("Known standardised to GNoME mix     %.2f%%   (difference %+.2f pp)"%(sk,cg-sk))
print("\nclass mix over these cells (GNoME / Known), %%:")
a=(j.groupby('anion_class').n_G.sum()/j.n_G.sum()*100)
b=(j.groupby('anion_class').n.sum()/j.n.sum()*100)
for k in a.index: print("   %-24s %5.1f  %5.1f"%(k,a[k],b[k]))
# bootstrap the standardised difference, resampling both sides independently
RNG=np.random.default_rng(20261011)
B=100000; d=[]
for _ in range(B):
    pg=RNG.binomial(j.n_G.values,j.flip_G.values)/j.n_G.values
    pk=RNG.binomial(j.n.values,j.flip_known.values)/j.n.values
    d.append(((pg-pk)*wK.values).sum()*100)
d=np.array(d)
print("\nstandardised GNoME-minus-Known difference at the Known mix: %+.2f pp [%.2f, %.2f]"
      %(d.mean(),np.percentile(d,2.5),np.percentile(d,97.5)))
# The LaTeX tables for Section 4.13 and Appendix D are written by
# phase37_matched_cluster.py, which uses the cluster bootstrap on the GNoME side. This
# script keeps the standardisation arithmetic and the all-binomial comparison only.
print("\n--- six further cells with 100 <= n_G < 500 (not in the Section 4.13 set) ---")
x=full[(full.n_G<500)&(full.n>=100)].sort_values("delta",ascending=False)
print(x[["n_elements","anion_class","has_Ln","has_3d","n_G","flip_G","n","flip_known","delta","lo","hi"]]
      .to_string(index=False,float_format=lambda v:f"{v:.3f}"))

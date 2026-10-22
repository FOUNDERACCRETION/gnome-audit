"""External review of 11 Oct 2026, priority 1: the matched-comparison intervals were an
independent-binomial parametric bootstrap on BOTH sides. On the GNoME side that assumes
independent materials and so carries none of the structural-prototype correction applied
everywhere else in the paper, while the paper claimed its rate intervals were
cluster-adjusted. With an intra-cluster correlation of 0.14 those intervals were too narrow.

This implements it properly:

  GNoME side  cluster bootstrap. One draw of prototypes with replacement per replicate,
              taken ACROSS all matched cells at once, so both within-cell dependence and
              dependence between cells sharing a prototype are preserved. Every cell rate
              and the standardised contrast are recomputed from that single draw.

  MP side     nonparametric material-level bootstrap. For a mean of 0/1 outcomes this is
              exactly binomial(n, p_hat)/n, which is what is drawn. No structural prototype
              is available for Materials Project entries -- the API fields we hold are
              material_id, energy_above_hull, nelements and chemsys -- so the MP side cannot
              be clustered and its interval remains an independence-based one. That is a
              stated limitation, not an oversight.

The two sides are drawn independently of each other because the datasets share no entries.
Reported for comparison: the old all-binomial interval, so the size of the correction is
visible rather than asserted.
"""
import ast, os, re
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
LN=set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
TM3D=set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
HAL=set("F Cl Br I".split()); CHAL=set("S Se Te".split()); PNIC=set("N P As Sb Bi".split())
def ac(e):
    if "O" in e: return "oxide"
    if e&HAL: return "halide"
    if e&CHAL: return "chalcogenide"
    if e&PNIC: return "pnictide"
    if "H" in e: return "hydride"
    if e&{"B","C","Si"}: return "boride/carbide/silicide"
    return "intermetallic"
def anon(c):
    cnt=[int(n) for _,n in re.findall(r"([A-Z][a-z]?)(\d+)",c)]
    return "".join(f"{chr(65+i)}{x}" for i,x in enumerate(sorted(cnt)))
NB=int(os.environ.get("NB","20000")); CH=int(os.environ.get("CH","500"))
SEED=int(os.environ.get("SEED","20261011"))
KEY=["n_elements","anion_class","has_Ln","has_3d"]

BASIS=os.environ.get("BASIS","release")   # "release" or "confirmed"
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["MaterialId","Composition","Elements","Space Group Number",C])
r[C]=pd.to_numeric(r[C],errors="coerce"); r=r.dropna(subset=[C]).reset_index(drop=True)
if BASIS=="confirmed":
    _s=pd.read_csv("data/stable_materials_summary.csv",usecols=["MaterialId",C])
    _s=_s.rename(columns={C:"pbe"}); _s["pbe"]=pd.to_numeric(_s.pbe,errors="coerce")
    r=r.merge(_s,on="MaterialId",how="left")
    r=r[r.pbe.notna()&(r.pbe<=0)].reset_index(drop=True)
print("basis: %s   GNoME rows %d"%(BASIS,len(r)),flush=True)
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_elements"]=els.map(len); r["anion_class"]=els.map(ac)
r["has_Ln"]=els.map(lambda e:bool(e&LN)); r["has_3d"]=els.map(lambda e:bool(e&TM3D))
r["flip"]=(r[C]>0).astype(float)
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)
kb=pd.read_csv("results/known_baseline.csv")
g=r.groupby(KEY).flip.agg(["size","mean"]).reset_index()
g.columns=KEY+["n_G","flip_G"]
j=kb.merge(g,on=KEY,how="inner")
j=j[(j.n_G>=500)&(j.n>=100)].copy()
j["delta"]=(j.flip_G-j.flip_known)*100
j=j.sort_values("delta",ascending=False).reset_index(drop=True)
NC=len(j)
print("matched cells: %d   GNoME entries %d   known entries %d"%(NC,j.n_G.sum(),j.n.sum()),flush=True)

# GNoME per-prototype, per-cell counts and flip sums, over the matched cells only
key_of={tuple(x):i for i,x in enumerate(j[KEY].itertuples(index=False,name=None))}
r["_cell"]=[key_of.get(t,-1) for t in zip(r.n_elements,r.anion_class,r.has_Ln,r.has_3d)]
sub=r[r._cell>=0]
codes,uniq=pd.factorize(sub.proto); P=len(uniq)
CNT=np.zeros((P,NC)); SUM=np.zeros((P,NC))
np.add.at(CNT,(codes,sub._cell.values),1.0)
np.add.at(SUM,(codes,sub._cell.values),sub.flip.values)
print("prototypes spanning the matched cells: %d   materials %d"%(P,len(sub)),flush=True)
spanning=(CNT>0).sum(axis=1)
print("prototypes appearing in more than one matched cell: %d (%.1f%%)"
      %((spanning>1).sum(),100*(spanning>1).mean()),flush=True)
nG=j.n_G.values.astype(float); nK=j.n.values.astype(float)
pG=j.flip_G.values; pK=j.flip_known.values
wK=nK/nK.sum(); wG=nG/nG.sum()

rng=np.random.default_rng(SEED)
pr=np.full(P,1.0/P)
DG=[]; DK=[]
done=0
while done<NB:
    k=min(CH,NB-done)
    W=rng.multinomial(P,pr,size=k).astype(float).T           # (P,k) prototype draw
    c=CNT.T@W; s=SUM.T@W                                     # (NC,k)
    with np.errstate(invalid="ignore",divide="ignore"):
        DG.append(np.where(c>0,s/np.where(c>0,c,1),np.nan))
    DK.append((rng.binomial(nK[:,None].astype(int).repeat(k,axis=1),
                            pK[:,None].repeat(k,axis=1))/nK[:,None]))
    done+=k
RG=np.concatenate(DG,axis=1); RK=np.concatenate(DK,axis=1)    # (NC, NB)
D=(RG-RK)*100
lo=np.nanpercentile(D,2.5,axis=1); hi=np.nanpercentile(D,97.5,axis=1)
# old all-binomial interval, for comparison only
rng2=np.random.default_rng(SEED)
BG=rng2.binomial(nG[:,None].astype(int).repeat(NB,axis=1),pG[:,None].repeat(NB,axis=1))/nG[:,None]
BK=rng2.binomial(nK[:,None].astype(int).repeat(NB,axis=1),pK[:,None].repeat(NB,axis=1))/nK[:,None]
Db=(BG-BK)*100
blo=np.percentile(Db,2.5,axis=1); bhi=np.percentile(Db,97.5,axis=1)
j["lo"]=lo; j["hi"]=hi; j["excl0"]=(lo>0)|(hi<0)
j["lo_binom"]=blo; j["hi_binom"]=bhi; j["excl0_binom"]=(blo>0)|(bhi<0)
j["widening"]=(hi-lo)/(bhi-blo)
print("\n%-34s %7s %8s %-20s %-20s %6s"%("cell","delta","","cluster 95% CI","all-binomial CI","widen"),flush=True)
NAME={"boride/carbide/silicide":"bor./carb./sil.","chalcogenide":"chalcogenide",
      "intermetallic":"intermetallic","oxide":"oxide","halide":"halide",
      "pnictide":"pnictide","hydride":"hydride"}
for _,x in j.iterrows():
    fl=[]
    if x.has_Ln: fl.append("Ln")
    if x.has_3d: fl.append("3d")
    nm="%d-el %s%s"%(x.n_elements,NAME[x.anion_class],(", "+", ".join(fl)) if fl else ", none")
    print("%-34s %+7.1f   %-20s %-20s %5.2f %s%s"
          %(nm,x.delta,"[%+.1f, %+.1f]"%(x.lo,x.hi),"[%+.1f, %+.1f]"%(x.lo_binom,x.hi_binom),
            x.widening,"*" if x.excl0 else " ","(was *)" if x.excl0_binom and not x.excl0 else ""),flush=True)
print("\nintervals excluding zero: cluster bootstrap %d of %d; all-binomial %d of %d"
      %(int(j.excl0.sum()),NC,int(j.excl0_binom.sum()),NC),flush=True)
print("median widening of the interval: %.2fx   (design effect implied: %.2f)"
      %(j.widening.median(),j.widening.median()**2),flush=True)
# standardised contrasts from the SAME prototype draw, weights held fixed
for lab,w in (("known-materials mix",wK),("GNoME mix",wG)):
    sd=np.nansum(D*w[:,None],axis=0)
    sb=np.sum(Db*w[:,None],axis=0)
    print("\nstandardised difference at the %-22s cluster %+.2f pp [%+.2f, %+.2f]   "
          "all-binomial %+.2f pp [%+.2f, %+.2f]"
          %(lab,np.nanmean(sd),np.nanpercentile(sd,2.5),np.nanpercentile(sd,97.5),
            sb.mean(),np.percentile(sb,2.5),np.percentile(sb,97.5)),flush=True)
j.to_csv("results/known_cells_clustered_%s.csv"%BASIS,index=False)
rows=[];sig=[]
for _,x in j.iterrows():
    fl=[]
    if x.has_Ln: fl.append("Ln")
    if x.has_3d: fl.append("3$d$")
    f=(", "+", ".join(fl)) if fl else ", no Ln/3$d$"
    line=("%d-el %s%s & %s & %.1f\\%% & %s & %.1f\\%% & $%+.1f$ & [$%+.1f$, $%+.1f$]"
          %(x.n_elements,NAME[x.anion_class],f,f"{int(x.n_G):,}".replace(",","{,}"),x.flip_G*100,
            f"{int(x.n):,}".replace(",","{,}"),x.flip_known*100,x.delta,x.lo,x.hi))
    rows.append(line+(" $\\ast$" if x.excl0 else "")+" \\\\")
    if x.excl0: sig.append(line+" \\\\")
if BASIS=="release":
    open("results/known_all_cells.tex","w").write("\n".join(rows)+"\n")
    open("results/known_tab11.tex","w").write("\n".join(sig)+"\n")
    print("\nwrote known_cells_clustered_release.csv, known_all_cells.tex (%d rows), known_tab11.tex (%d rows)"
          %(len(rows),len(sig)),flush=True)
else:
    print("\nwrote known_cells_clustered_confirmed.csv",flush=True)

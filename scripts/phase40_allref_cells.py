"""Headline composition cells under the undocumented "All" reference set.

Section 4.16 reports the reference-set sensitivity at anion-class level only, which answers
whether the broad ordering survives but not whether the headline cell contrast does. This
recomputes the two primary cells -- lanthanide-free 3d-containing quaternary oxides and
lanthanide-free 3d-free quaternary intermetallics -- against the "All" column, with the
spread and the ratio drawn jointly so that prototypes spanning both cells are handled.

Cell definitions, prototype key and the cluster bootstrap are identical to
phase34_converged_intervals.py; only the outcome column changes.
"""
import ast, os, re
import numpy as np, pandas as pd
C="Decomposition Energy Per Atom MP"
CA="Decomposition Energy Per Atom All"
CRel="Decomposition Energy Per Atom Relative"
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

NB=int(os.environ.get("NB","20000")); CH=int(os.environ.get("CH","1000"))
SEED=int(os.environ.get("SEED","20261009"))
r=pd.read_csv("data/stable_materials_r2scan.csv",
              usecols=["Composition","Elements","Space Group Number",C,CA,CRel])
for col in (C,CA,CRel): r[col]=pd.to_numeric(r[col],errors="coerce")
els=r["Elements"].map(ast.literal_eval).map(set)
r["n_el"]=els.map(len); r["anion"]=els.map(ac)
r["Ln"]=els.map(lambda e:bool(e&LN)); r["d3"]=els.map(lambda e:bool(e&TM3D))
r["proto"]=r["Space Group Number"].astype(str)+"|"+r.Composition.map(anon)

print(f"rows {len(r)}; non-null MP {r[C].notna().sum()}, non-null All {r[CA].notna().sum()}")
a=r.dropna(subset=[CA]).reset_index(drop=True)
print(f"global All flip rate {100*(a[CA]>0).mean():.1f}%  (paper: 25.8%)")
# The release documents a "Relative" column as the distance to the hull of all entries
# except the current one -- the leave-one-out phase-separation energy. It is EMPTY in both
# released files, so it cannot be used as an outcome. Asserted here so a future run notices
# if DeepMind ever populates it.
nrel=int(pd.to_numeric(r[CRel],errors="coerce").notna().sum())
print(f"Relative column populated rows: {nrel} of {len(r)} (documented but empty)")
assert nrel==0, "the Relative column is now populated -- reconsider it as the primary outcome"

CELLS={"oxide, no Ln, 3d":      dict(n_el=4,anion="oxide",        Ln=False,d3=True),
       "intermetallic, no Ln/3d":dict(n_el=4,anion="intermetallic",Ln=False,d3=False)}

def mask(d,spec):
    m=np.ones(len(d),bool)
    for k,v in spec.items(): m&=(d[k]==v).values
    return m

def joint(d,col,specs,seed=SEED,B=NB,ch=CH):
    """One prototype draw per replicate, shared across every cell in specs."""
    flip=(d[col]>0).astype(float).values
    codes,uniq=pd.factorize(d.proto); G=len(uniq)
    N=np.zeros((G,len(specs))); S=np.zeros((G,len(specs)))
    for j,spec in enumerate(specs):
        m=mask(d,spec)
        np.add.at(N[:,j],codes[m],1.0)
        np.add.at(S[:,j],codes[m],flip[m])
    rng=np.random.default_rng(seed); p=np.full(G,1.0/G); out=[]; done=0
    while done<B:
        k=min(ch,B-done)
        W=rng.multinomial(G,p,size=k).astype(float).T          # (G,k)
        c=N.T@W; s=S.T@W                                        # (cells,k)
        out.append(np.where(c>0,s/np.where(c>0,c,1.0),np.nan))
        done+=k
    V=np.concatenate(out,axis=1)*100.0                          # (cells,B) in %
    pts=[100.0*S[:,j].sum()/N[:,j].sum() for j in range(len(specs))]
    ns=[int(N[:,j].sum()) for j in range(len(specs))]
    gs=[int((N[:,j]>0).sum()) for j in range(len(specs))]
    return np.array(pts),ns,gs,V

names=list(CELLS); specs=[CELLS[k] for k in names]
print(f"\n{NB} replicates, seed {SEED}, cluster bootstrap over structural prototypes\n")
rowsout=[]
for col,lab in ((C,"MP"),(CA,"All")):
    d=r.dropna(subset=[col]).reset_index(drop=True)
    pts,ns,gs,V=joint(d,col,specs)
    lo=np.percentile(V,2.5,axis=1); hi=np.percentile(V,97.5,axis=1)
    sp=V[0]-V[1]; ra=V[0]/V[1]
    print(f"--- {lab} reference ---")
    for j,nm in enumerate(names):
        print(f"  {nm:26s} n={ns[j]:6d} clusters={gs[j]:5d}  "
              f"{pts[j]:5.1f}%  [{lo[j]:.1f}, {hi[j]:.1f}]")
        rowsout.append(dict(reference=lab,cell=nm,n=ns[j],n_clusters=gs[j],
                            rate_pct=round(pts[j],1),lo=round(lo[j],1),hi=round(hi[j],1)))
    s_pt=pts[0]-pts[1]; r_pt=pts[0]/pts[1]
    print(f"  spread {s_pt:+.1f} pp  [{np.percentile(sp,2.5):+.1f}, {np.percentile(sp,97.5):+.1f}]")
    print(f"  ratio  {r_pt:.2f}x     [{np.percentile(ra,2.5):.2f}, {np.percentile(ra,97.5):.2f}]")
    rowsout.append(dict(reference=lab,cell="spread (oxide - intermetallic, pp)",n="",n_clusters="",
                        rate_pct=round(s_pt,1),lo=round(np.percentile(sp,2.5),1),
                        hi=round(np.percentile(sp,97.5),1)))
    rowsout.append(dict(reference=lab,cell="ratio (oxide / intermetallic)",n="",n_clusters="",
                        rate_pct=round(r_pt,2),lo=round(np.percentile(ra,2.5),2),
                        hi=round(np.percentile(ra,97.5),2)))
pd.DataFrame(rowsout).to_csv("results/allref_cells.csv",index=False)
print("\nwrote results/allref_cells.csv")

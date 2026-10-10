"""Two referee-flagged analyses on the Materials Project side, stdlib + requests only.
(1) Hull-completeness probe: for each GGA-stable MP material, what fraction of the GGA
    entries in its chemical system and all sub-systems have an r2SCAN calculation?
    Then: does the chemistry flip gradient survive within completeness strata?
(2) Hull-distance-matched GNoME-vs-MP comparison: compare flip rates within
    (anion, 3d) x PBE-hull-distance bins, not chemistry alone.
Run from the repo root:  python3 scripts/q2b_completeness.py"""
import csv, json, time, random, ast
from collections import defaultdict
from itertools import combinations
import requests
random.seed(20261009)
import os
KEY=os.environ.get("MP_API_KEY") or open(os.environ.get("MP_API_KEY_FILE","/mnt/c/research/mp_api_key.txt")).read().strip()
H={"X-API-KEY":KEY}; B="https://api.materialsproject.org/materials/thermo/"
OUT="results"
FIELDS="material_id,thermo_type,energy_above_hull,nelements,chemsys"
def pull(tt):
    rows=[]; page=1
    while True:
        r=None
        for a in range(4):
            try:
                r=requests.get(B,params={"thermo_types":tt,"_per_page":1000,"_page":page,
                   "_fields":FIELDS},headers=H,timeout=120)
                if r.status_code==200: break
                time.sleep(3)
            except Exception as e: time.sleep(5); r=None
        if r is None or r.status_code!=200: break
        d=r.json().get("data",[])
        if not d: break
        rows+=d
        if page%40==0: print(f"  {tt} {len(rows)}",flush=True)
        page+=1
    return rows
print("pulling r2SCAN",flush=True); sc=pull("r2SCAN")
print("pulling GGA_GGA+U",flush=True); gg=pull("GGA_GGA+U")
with open(f"{OUT}/mp_r2scan_raw.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["material_id","energy_above_hull","nelements","chemsys"]);w.writeheader()
    for d in sc: w.writerow({k:d.get(k) for k in w.fieldnames})
with open(f"{OUT}/mp_gga_raw.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["material_id","energy_above_hull","nelements","chemsys"]);w.writeheader()
    for d in gg: w.writerow({k:d.get(k) for k in w.fieldnames})
print(f"saved raw: r2SCAN {len(sc)}, GGA {len(gg)}",flush=True)

scan={d["material_id"]:d["energy_above_hull"] for d in sc if d.get("energy_above_hull") is not None}
# index GGA entries by chemsys
by_cs=defaultdict(list)
for d in gg:
    if d.get("chemsys"): by_cs[d["chemsys"]].append(d["material_id"])
def subsystems(cs):
    el=cs.split("-"); out=[]
    for k in range(1,len(el)+1):
        for c in combinations(el,k): out.append("-".join(sorted(c)))
    return out
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
TOL=1e-6
rows=[]
for d in gg:
    mid=d["material_id"]; eg=d.get("energy_above_hull"); cs=d.get("chemsys")
    if eg is None or mid not in scan or not cs or eg>TOL: continue
    comp=[m for s in subsystems(cs) for m in by_cs.get(s,[]) if m!=mid]
    cov=sum(1 for m in comp if m in scan)/len(comp) if comp else float("nan")
    e=set(cs.split("-"))
    rows.append(dict(material_id=mid,chemsys=cs,n_el=len(e),anion=ac(e),
        d3=bool(e&TM3D),Ln=bool(e&LN),n_competitors=len(comp),completeness=cov,
        flip=1 if scan[mid]>TOL else 0))
with open(f"{OUT}/mp_completeness.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"\nGGA-stable MP materials with r2SCAN: {len(rows)}",flush=True)

# (1a) completeness by anion class
print("\n=== (1a) r2SCAN hull completeness by anion class (MP side) ===")
byA=defaultdict(list)
for r in rows:
    if r["completeness"]==r["completeness"]: byA[r["anion"]].append(r)
print(f"{'anion':24s}{'n':>7}{'median compl.':>15}{'flip':>8}")
for a,v in sorted(byA.items(),key=lambda kv:-sum(x['flip'] for x in kv[1])/len(kv[1])):
    c=sorted(x["completeness"] for x in v); med=c[len(c)//2]
    print(f"{a:24s}{len(v):7d}{med*100:14.1f}%{sum(x['flip'] for x in v)/len(v)*100:7.1f}%")

# (1b) chemistry gradient WITHIN completeness strata
print("\n=== (1b) flip rate by completeness tercile x (3d present) ===")
cs_sorted=sorted(r["completeness"] for r in rows if r["completeness"]==r["completeness"])
t1=cs_sorted[len(cs_sorted)//3]; t2=cs_sorted[2*len(cs_sorted)//3]
def terc(c): return "low" if c<t1 else ("mid" if c<t2 else "high")
cell=defaultdict(list)
for r in rows:
    if r["completeness"]!=r["completeness"]: continue
    cell[(terc(r["completeness"]),r["d3"])].append(r["flip"])
print(f"terciles cut at completeness {t1*100:.1f}% and {t2*100:.1f}%")
print(f"{'tercile':8s}{'3d':>6}{'n':>7}{'flip':>8}")
for k in sorted(cell): v=cell[k]; print(f"{k[0]:8s}{str(k[1]):>6}{len(v):7d}{sum(v)/len(v)*100:7.1f}%")
# and the oxide+3d vs intermetallic contrast within each tercile
cell2=defaultdict(list)
for r in rows:
    if r["completeness"]!=r["completeness"]: continue
    lab=("oxide+3d" if r["anion"]=="oxide" and r["d3"] else
         "interm" if r["anion"]=="intermetallic" and not r["d3"] and not r["Ln"] else None)
    if lab: cell2[(terc(r["completeness"]),lab)].append(r["flip"])
print("\n   oxide+3d vs plain intermetallic within each tercile:")
for tc in ("low","mid","high"):
    a=cell2.get((tc,"oxide+3d"),[]); b=cell2.get((tc,"interm"),[])
    if len(a)>=30 and len(b)>=30:
        print(f"   {tc:5s}  oxide+3d {sum(a)/len(a)*100:5.1f}% (n={len(a)})   "
              f"interm {sum(b)/len(b)*100:5.1f}% (n={len(b)})   gap {(sum(a)/len(a)-sum(b)/len(b))*100:+5.1f}pp")

# (2) hull-distance-matched GNoME vs MP
print("\n=== (2) GNoME vs MP, matched on (anion,3d) AND PBE hull-distance bin ===")
C="Decomposition Energy Per Atom MP"
gs={}
with open("data/stable_materials_summary.csv") as f:
    for r_ in csv.DictReader(f):
        try: gs[r_["MaterialId"]]=float(r_[C])
        except: pass
gcell=defaultdict(list)
with open("data/stable_materials_r2scan.csv") as f:
    for r_ in csv.DictReader(f):
        mid=r_["MaterialId"]
        if mid not in gs or gs[mid]>0: continue
        try: sv=float(r_[C])
        except: continue
        e=set(ast.literal_eval(r_["Elements"]))
        if len(e)!=3: continue      # ternaries: where MP cells are largest
        pbe=-gs[mid]                 # distance below hull, positive
        b="near" if pbe<0.02 else ("mid" if pbe<0.05 else "far")
        gcell[(ac(e),bool(e&TM3D),b)].append(1 if sv>0 else 0)
# MP side needs GGA distance-below-hull — but GGA-stable means e_hull==0 exactly.
# MP has no 'how far below hull' for stable entries (it's a hull energy, floored at 0).
print("   NOTE: MP energy_above_hull is 0 for every stable entry; MP exposes no 'distance")
print("   below hull' (decomposition energy) in this endpoint. Matching on hull distance")
print("   therefore cannot be done symmetrically. Reporting GNoME-side stratification only:")
print(f"   {'cell':30s}{'near(<20meV)':>14}{'mid':>10}{'far(>50)':>10}")
seen=set()
for (a,d3,b) in gcell:
    if (a,d3) in seen: continue
    seen.add((a,d3))
    vals=[]
    for bb in ("near","mid","far"):
        v=gcell.get((a,d3,bb),[])
        vals.append(f"{sum(v)/len(v)*100:5.1f}% ({len(v)})" if len(v)>=100 else "     ---")
    print(f"   {a[:16]+('+3d' if d3 else ''):30s}{vals[0]:>14}{vals[1]:>10}{vals[2]:>10}")
print("\nDONE")

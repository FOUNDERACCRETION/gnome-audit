"""Q2 known-materials baseline. Standard library only (plus requests)."""
import csv, json, random, time
from collections import defaultdict
import requests

import os
KEY=os.environ.get("MP_API_KEY") or open(os.environ.get("MP_API_KEY_FILE","/mnt/c/research/mp_api_key.txt")).read().strip()
H={"X-API-KEY":KEY}; B="https://api.materialsproject.org/materials/thermo/"
OUT="results"
FIELDS="material_id,thermo_type,energy_above_hull,nelements,chemsys"
random.seed(20261009)

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
            except Exception as e:
                print(f"  retry {a}: {e}",flush=True); time.sleep(5); r=None
        if r is None or r.status_code!=200:
            print(f"  giving up at page {page}",flush=True); break
        d=r.json(); data=d.get("data",[])
        if not data: break
        rows+=data
        if page==1 or page%25==0:
            print(f"  {tt} page {page}: {len(rows)} rows of {d.get('meta',{}).get('total_doc')}",flush=True)
        page+=1
    return rows

print("pulling r2SCAN...",flush=True); sc=pull("r2SCAN"); print(f"r2SCAN {len(sc)}",flush=True)
print("pulling GGA_GGA+U...",flush=True); gg=pull("GGA_GGA+U"); print(f"GGA {len(gg)}",flush=True)

scan={d["material_id"]:d["energy_above_hull"] for d in sc if d.get("energy_above_hull") is not None}
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
cells=defaultdict(list); matched=0; tot_flip=0
for d in gg:
    mid=d["material_id"]; eg=d.get("energy_above_hull")
    if eg is None or mid not in scan: continue
    matched+=1
    if eg>TOL: continue                      # only GGA-stable materials can flip
    es=scan[mid]; flip=1 if es>TOL else 0; tot_flip+=flip
    e=set(str(d.get("chemsys","")).split("-"))
    cells[(d.get("nelements"),ac(e),bool(e&LN),bool(e&TM3D))].append(flip)

n_stable=sum(len(v) for v in cells.values())
print(f"\nmatched both functionals: {matched}",flush=True)
print(f"GGA-stable subset: {n_stable}  flip rate "
      f"{(tot_flip/n_stable*100) if n_stable else 0:.2f}%",flush=True)

def ci(v,nb=1000):
    n=len(v); p=sum(v)/n
    d=sorted(sum(1 for _ in range(n) if random.random()<p)/n for _ in range(nb))
    return p,d[int(.025*nb)],d[int(.975*nb)]

rows=[]
for (ne,a,ln,d3),v in cells.items():
    if len(v)<100: continue
    p,lo,hi=ci(v)
    rows.append(dict(n_elements=ne,anion_class=a,has_Ln=ln,has_3d=d3,
                     n=len(v),flip_known=round(p,4),lo=round(lo,4),hi=round(hi,4)))
rows.sort(key=lambda r:-r["flip_known"])
with open(f"{OUT}/known_baseline.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\n=== KNOWN MATERIALS: flip by cell (n>=100) ===",flush=True)
print(f"{'n_el':>4} {'anion':<24} {'Ln':>5} {'3d':>5} {'n':>7} {'flip':>7} {'lo':>7} {'hi':>7}",flush=True)
for r in rows:
    print(f"{r['n_elements']:>4} {r['anion_class']:<24} {str(r['has_Ln']):>5} "
          f"{str(r['has_3d']):>5} {r['n']:>7} {r['flip_known']:>7.4f} {r['lo']:>7.4f} {r['hi']:>7.4f}",flush=True)
print("\nwrote known_baseline.csv",flush=True)

import json, requests
KEY=open("/mnt/c/research/mp_api_key.txt").read().strip()
H={"X-API-KEY":KEY}; B="https://api.materialsproject.org/materials/thermo/"
def go(label, params):
    params.setdefault("_per_page",3); params.setdefault("_page",1)
    params.setdefault("_fields","material_id,thermo_type,energy_above_hull,formula_pretty")
    try:
        r=requests.get(B,params=params,headers=H,timeout=60)
        m=r.json().get("meta",{}) if r.status_code==200 else {}
        print(f"{label:48s} HTTP {r.status_code}  total_doc={m.get('total_doc')}")
        if r.status_code==200:
            for d in r.json().get("data",[])[:2]: print("     ",d)
        else: print("     ",r.text[:200])
    except Exception as e: print(f"{label:48s} FAIL {e}")

go("thermo_types=R2SCAN",            {"thermo_types":"R2SCAN"})
go("thermo_types=GGA_GGA+U_R2SCAN",  {"thermo_types":"GGA_GGA+U_R2SCAN"})
go("chemsys=Pr + R2SCAN",            {"chemsys":"Pr","thermo_types":"R2SCAN"})
go("chemsys=Pr (no filter)",         {"chemsys":"Pr"})
go("no thermo_types at all",         {})
print("\n--- what thermo_type values actually appear in a sample of 200 docs:")
r=requests.get(B,params={"_per_page":200,"_page":1,"_fields":"thermo_type"},headers=H,timeout=60)
from collections import Counter
print(Counter(d["thermo_type"] for d in r.json().get("data",[])))

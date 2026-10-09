"""10-second smoke test: can we reach the MP thermo endpoint and what comes back?"""
import json, requests, sys
KEY = open("/mnt/c/research/mp_api_key.txt").read().strip()
H = {"X-API-KEY": KEY}
base = "https://api.materialsproject.org/materials/thermo/"
for tt in ("GGA_GGA+U", "R2SCAN"):
    p = {"thermo_types": tt, "_per_page": 5, "_page": 1,
         "_fields": "material_id,thermo_type,energy_above_hull,formula_pretty,nelements"}
    try:
        r = requests.get(base, params=p, headers=H, timeout=60)
        print(f"\n--- {tt}: HTTP {r.status_code}")
        if r.status_code == 200:
            d = r.json()
            print("meta:", json.dumps(d.get("meta", {}))[:300])
            for row in d.get("data", [])[:3]:
                print("  ", row)
        else:
            print(r.text[:400])
    except Exception as e:
        print(f"{tt}: FAILED {type(e).__name__}: {e}")

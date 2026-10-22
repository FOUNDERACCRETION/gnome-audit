"""Pull the Materials Project sub-hull margin, so the matched comparison can be depth-matched.

WHY THIS IS A SEPARATE SCRIPT. `api.materialsproject.org` is unreachable from the analysis
container and from the Linux VM used for this work: both get HTTP 403 at the egress proxy
on CONNECT. This script therefore has to run on the machine that holds the API key and has
outbound access -- the same machine that ran q2_known_baseline.py.

WHAT IT ADDS. The earlier pull took four fields, of which `energy_above_hull` is floored at
zero for stable entries and so carries no information about how far BELOW the hull a stable
entry sits. Two further fields do:

  decomposition_enthalpy                 the phase-separation energy of Bartel et al.,
                                         obtained by removing the entry from its own phase
                                         diagram. Negative for stable entries. Populated
                                         for all entries. May be None if the optimiser
                                         failed, so handle None.
  equilibrium_reaction_energy_per_atom   the same margin, set only where is_stable is True.

RUN IT:
    python scripts/phase42_pull_hull_depth.py

It reads the key from mp_api_key.txt (never printed), writes
results/mp_gga_hulldepth.csv and results/mp_r2scan_hulldepth.csv, and prints a sanity
check that the margin is negative for stable entries. Nothing else in the pipeline is
touched until phase43 consumes these files.
"""
import os, sys
from pathlib import Path

KEYFILE = Path("mp_api_key.txt")
if not KEYFILE.exists():
    sys.exit("mp_api_key.txt not found; run this from the repository root on the machine "
             "that holds the key")
KEY = KEYFILE.read_text().strip()          # deliberately never printed or logged

try:
    from mp_api.client import MPRester
except ImportError:
    sys.exit("pip install mp-api")

# The stored value for the pure r2SCAN hull is spelled "r2SCAN". The ThermoType enum spells
# it "R2SCAN" and recent clients coerce it; pass the enum where available so either works.
FIELDS = ["material_id", "energy_above_hull", "decomposition_enthalpy",
          "equilibrium_reaction_energy_per_atom", "is_stable", "nelements", "chemsys"]

import pandas as pd
with MPRester(KEY) as mpr:
    for tt, out in (("GGA_GGA+U", "results/mp_gga_hulldepth.csv"),
                    ("R2SCAN",    "results/mp_r2scan_hulldepth.csv")):
        print(f"querying thermo_types={tt} ...", flush=True)
        docs = mpr.materials.thermo.search(thermo_types=[tt], fields=FIELDS)
        rows = [{f: getattr(d, f, None) for f in FIELDS} for d in docs]
        df = pd.DataFrame(rows)
        df.to_csv(out, index=False)
        st = df[df.is_stable == True]
        print(f"  {len(df)} docs -> {out}")
        print(f"  stable entries: {len(st)}")
        for c in ("decomposition_enthalpy", "equilibrium_reaction_energy_per_atom"):
            v = pd.to_numeric(st[c], errors="coerce").dropna()
            if len(v):
                print(f"  {c}: {len(v)} populated, "
                      f"{100*(v <= 0).mean():.1f}% non-positive, median {v.median():.4f}")
            else:
                print(f"  {c}: EMPTY for stable entries -- do not use it as the covariate")
print("\nSANITY CHECK: the covariate is usable only if one of the two fields above is "
      "populated and non-positive for essentially every stable entry.")

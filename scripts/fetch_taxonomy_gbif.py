"""
Query GBIF taxonomy for all species in zims_all.csv.
Writes results/zims_taxonomy.csv with genus, family, order columns.
Run once; the combined viewer build script reads this file.
"""
import requests, pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

ZIMS_CSV = "results/zims_all.csv"
OUT_CSV  = "results/zims_taxonomy.csv"

df = pd.read_csv(ZIMS_CSV)
species_list = list(df["binSpecies"].unique())
print(f"Querying GBIF for {len(species_list)} species...")

def fetch(species):
    try:
        r = requests.get("https://api.gbif.org/v1/species/match",
                         params={"name": species, "verbose": False}, timeout=10)
        d = r.json()
        return species, {"genus": d.get("genus",""), "family": d.get("family",""),
                         "order": d.get("order",""), "confidence": d.get("confidence",0)}
    except:
        return species, {"genus":"","family":"","order":"","confidence":0}

results = {}
with ThreadPoolExecutor(max_workers=12) as ex:
    futures = {ex.submit(fetch, sp): sp for sp in species_list}
    for i, fut in enumerate(as_completed(futures)):
        sp, tax = fut.result()
        results[sp] = tax
        if (i+1) % 100 == 0:
            print(f"  {i+1}/{len(species_list)}")

rows = [{"binSpecies": sp, **tax} for sp, tax in results.items()]
pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f"Saved {OUT_CSV}  ({len(rows)} species, "
      f"{pd.DataFrame(rows)['family'].nunique()} families, "
      f"{pd.DataFrame(rows)['order'].nunique()} orders)")

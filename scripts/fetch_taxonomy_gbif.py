"""
Query GBIF taxonomy for all species in zims_all.csv.
Writes results/zims_taxonomy.csv with genus, family, order, common_name columns.
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
        key = d.get("usageKey", "")
        common = ""
        if key:
            vr = requests.get(f"https://api.gbif.org/v1/species/{key}/vernacularNames",
                              params={"limit": 20}, timeout=10)
            names = [x["vernacularName"] for x in vr.json().get("results", [])
                     if x.get("language", "") == "eng"]
            if names:
                common = names[0]
        return species, {
            "genus": d.get("genus", ""),
            "family": d.get("family", ""),
            "order": d.get("order", ""),
            "confidence": d.get("confidence", 0),
            "common_name": common,
        }
    except Exception:
        return species, {"genus": "", "family": "", "order": "",
                         "confidence": 0, "common_name": ""}

results = {}
with ThreadPoolExecutor(max_workers=16) as ex:
    futures = {ex.submit(fetch, sp): sp for sp in species_list}
    for i, fut in enumerate(as_completed(futures)):
        sp, tax = fut.result()
        results[sp] = tax
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(species_list)}")

rows = [{"binSpecies": sp, **tax} for sp, tax in results.items()]
out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)
n_with = (out["common_name"] != "").sum()
print(f"Saved {OUT_CSV}  ({len(rows)} species, {n_with} with English common names, "
      f"{out['family'].nunique()} families, {out['order'].nunique()} orders)")

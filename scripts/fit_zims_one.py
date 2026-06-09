"""
Fit SR FP model to one (class, species, sex) from the ZIMS dataset.
- Xc=1, k free (κ free), ndims=4; ndims=5 (Makeham external hazard) fallback if rms > RMS_THRESH
- Infant deaths excluded: lifespan < 0.10 * L_prelim (estimated from dead animals)
Usage: python fit_zims_one.py <class> <binSpecies> <sex_code>
  sex_code: f or m
"""
import sys, os, math
import numpy as np
import pandas as pd

WORK      = "/home/labs/alon/avimayo/sr_fits"
DATA      = f"{WORK}/data/zims/zims_clean.csv"
OUTDIR    = f"{WORK}/results/zims"
RMS_THRESH = 0.05   # if ndims=4 rms exceeds this, also try ndims=5

sys.path.insert(0, f"{WORK}/scripts")
sys.path.insert(0, f"{WORK}/scripts/kfix")
from fp_curvefit_seed    import fp_curvefit_seed
from fp_curvefit_makeham import fp_curvefit_makeham

cls      = sys.argv[1]
species  = sys.argv[2]
sex_code = sys.argv[3]   # "f" or "m"

sex_label  = "Female" if sex_code == "f" else "Male"
safe_sp    = species.replace(" ", "_").replace("/", "-")
outfile    = f"{OUTDIR}/{cls}_{safe_sp}_{sex_code}.csv"

if os.path.exists(outfile):
    print(f"Already done: {outfile}")
    sys.exit(0)

df  = pd.read_csv(DATA)
sub = df[(df["binSpecies"] == species) & (df["SexType"] == sex_label)].copy()

# ── infant-death cutoff: 10% of preliminary mean lifespan ────────────────────
dead_times   = sub.loc[sub["event"] == 1, "lifespan_days"].values
L_prelim     = dead_times.mean() if len(dead_times) > 0 else sub["lifespan_days"].mean()
infant_cutoff = 0.10 * L_prelim
sub = sub[sub["lifespan_days"] >= infant_cutoff].copy()

death_times = sub["lifespan_days"].values.astype(float)
events      = sub["event"].values.astype(int)

n      = len(death_times)
n_dead = int(events.sum())
print(f"{cls} | {species} | {sex_code}  n={n}  n_dead={n_dead}  infant_cutoff={infant_cutoff:.1f}d", flush=True)

# ── ndims=4 fit ───────────────────────────────────────────────────────────────
r4 = fp_curvefit_seed(death_times, events=events, xc=1.0, ndims=4)

result = r4
ndims_used = 4

# ── ndims=5 fallback if rms is poor ──────────────────────────────────────────
if r4["residual"] > RMS_THRESH:
    try:
        r5 = fp_curvefit_makeham(death_times, events=events, xc=1.0)
        if r5["residual"] < r4["residual"]:
            result = r5
            ndims_used = 5
            print(f"  → ndims=5 better: rms {r4['residual']:.4f} → {r5['residual']:.4f}", flush=True)
        else:
            print(f"  → ndims=4 kept: rms4={r4['residual']:.4f} rms5={r5['residual']:.4f}", flush=True)
    except Exception as e:
        print(f"  → ndims=5 failed ({e}), keeping ndims=4", flush=True)

s           = math.sqrt(result["rho_eta"]) / result["rho_eps"]
beta_xc_eps = result["rho_beta"] / result["rho_eps"]

row = {
    "class":            cls,
    "binSpecies":       species,
    "sex":              sex_code,
    "n":                n,
    "n_dead":           n_dead,
    "infant_cutoff_d":  infant_cutoff,
    "ndims":            ndims_used,
    "L":                result["L"],
    "rms":              result["residual"],
    "arm":              result["arm"],
    "wall_s":           result["wall_s"],
    "rho_eta":          result["rho_eta"],
    "rho_beta":         result["rho_beta"],
    "rho_eps":          result["rho_eps"],
    "kappa":            result["kappa"],
    "omega":            result["omega"],
    "s":                s,
    "beta_xc_eps":      beta_xc_eps,
    "external_hazard":  result.get("external_hazard", float("nan")),
}

os.makedirs(OUTDIR, exist_ok=True)
pd.DataFrame([row]).to_csv(outfile, index=False)
print(f"Saved {outfile}", flush=True)

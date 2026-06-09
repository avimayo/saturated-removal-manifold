# Saturated-removal lifetime manifold — interactive viewer

An interactive 3-D view of the lifetime **constraint surface** of the saturated-removal (SR) aging model,
with fitted species placed in the dimensionless coordinates

$$\rho_\beta = L\beta/X_c, \qquad \rho_\eta = L^2\eta/X_c, \qquad \rho_\epsilon = L\epsilon/X_c^2,$$

formed from the mean lifetime $L$, the lethal threshold $X_c$, the maximal removal rate $\beta$,
the production-acceleration $\eta$, and the diffusion strength $\epsilon$.

## ▶ Open the viewers

| File | Description |
|------|-------------|
| [`index.html`](https://avimayo.github.io/saturated-removal-manifold/) | Original manifold — fitted multi-species points on the surface |
| [`viewer.html`](https://raw.githack.com/avimayo/saturated-removal-manifold/main/viewer.html) | **Combined viewer** — manifold + ITP interventions + 1 633 zoo-animal fits + 2-D analysis + phylo tab |

> **Note** `viewer.html` is ~15 MB; download and open locally if the CDN link is slow.

## What you're looking at

- The tan **surface** is the mean-lifetime constraint — a two-dimensional manifold in the three rates.
  Requiring the model's mean first-passage lifetime to equal the observed one gives one closed-form relation,
  $\tfrac12\rho_\eta-\rho_\beta = 1 - 1/\omega + e^{-\omega}/\omega$ with $\omega=(\tfrac12\rho_\eta-\rho_\beta)/\rho_\epsilon$.
- The dark **ridge** is the one-dimensional curve the surface singles out in closed form,
  labelled by $\omega$ (increasing toward the production end).
- **Points** are SR-fitted species/interventions; they fall along the ridge.

## Data layers in `viewer.html`

### NIA Interventions Testing Program (ITP) — mouse interventions

FP-fitted survival curves from the NIA ITP: 64 curves across 28 interventions (rapamycin, 17α-E2,
metformin, and others), 5 cohorts (C2009–C2016), male and female.
The fit is κ-free with X_c = 1.

**Results:** [`results/itp_fp_fits.csv`](results/itp_fp_fits.csv)

### ZIMS zoo-animal survey — 1 633 species × sex curves

SR FP-fitter applied to the ZIMS (Species360) zoo-animal dataset: 1.4 M animals across
Mammalia, Aves, Reptilia, Amphibia, and Chondrichthyes.

**Inclusion criteria:**
- Exact known birth date (`BirthDateEstimateType` blank) — excludes 173 686 animals with
  estimated or range birth dates.
- Lifespan ≥ 1 day — excludes 30 399 animals with a departure/death date on or before birth.
- Infant deaths excluded per curve: lifespan < 10 % of the species×sex mean dead lifespan.
- ≥ 50 confirmed deaths per species × sex to attempt a fit.

**Fit:** κ free, X_c = 1, ndims = 4; ndims = 5 (Makeham external hazard) used as a fallback
when ndims = 4 RMS > 0.05.  87 / 1 633 curves (5.3 %) needed ndims = 5.

**Results:** [`results/zims_all.csv`](results/zims_all.csv)  
**Taxonomy:** [`results/zims_taxonomy.csv`](results/zims_taxonomy.csv) (GBIF — genus / family / order)

### Goodness of fit summary

| Metric | Value |
|--------|-------|
| Median RMS | 0.024 |
| Curves with RMS > 0.10 | 3 (*Cacatua leadbeateri* ♂, *Cacatua moluccensis* ♂, *Ovis aries* ♂) |
| Arm: removal / production | 95 % / 5 % |
| ndims = 5 used | 87 (5.3 %) |

## Combined viewer — tabs and filters

**Tab 1 — 3D Manifold**  
Manifold surface + ITP points + ZIMS points. Right-side filter panel:
- *ITP*: intervention dropdown (with category-level selection) + sex toggle.
- *ZIMS*: class checkboxes with emoji (🦁🦅🦎🐸🦈), sex, arm, RMS quality slider, genus dropdown.

**Tab 2 — 2D Analysis**  
Interactive scatter panels: sharpness *s* vs mean lifespan *L* (with ITP ★ overlaid),
and identifiable parameter space β·X_c/ε vs *s*, both coloured by taxonomic class.

**Tab 3 — Phylo Analysis**  
Violin distributions of Euclidean distance in log(ρ_β, ρ_η, ρ_ε)-space grouped by
taxonomic proximity (same genus → family → order → class → cross-class).
Same-genus species are closer in parameter space (median distance 0.94)
than cross-class pairs (1.24), but the overlap between levels is large.

> **Manifold caveat:** the surface was computed at κ = 0.
> Naveh's original fits also use κ = 0. The ITP and ZIMS FP fits use κ free (X_c = 1),
> so those points need not lie exactly on the surface.

## Reproducing the ZIMS fits

Fits were run on the WEXAC cluster (Weizmann Institute, LSF scheduler):

```bash
# 1. Prepare data locally → results/zims_all.csv
#    (requires the ZIMS per-animal CSV files, not included in this repo)

# 2. Upload to WEXAC and submit array job
scp results/zims_all.csv     wexac:~/sr_fits/data/zims/
scp scripts/fit_zims_one.py  wexac:~/sr_fits/scripts/zims/
scp scripts/run_zims_array.sh wexac:~/sr_fits/scripts/zims/

N=$(wc -l < ~/sr_fits/zims_groups.txt)
bsub -q short -J "zims_fp[1-$N]" -n 1 \
  -R "rusage[mem=2000] span[hosts=1]" \
  -o ~/sr_fits/logs/zims/%J_%I.out \
  bash ~/sr_fits/scripts/zims/run_zims_array.sh
```

Python environment: `/home/labs/alon/navehr/.conda/envs/srtools/bin/python`  
(requires `fp_curvefit_seed.py`, `fp_curvefit_makeham.py`, `fpt_full_model.py` from the FP-fitter repo)

## Rebuilding `viewer.html`

```bash
# Optionally refresh taxonomy
python scripts/fetch_taxonomy_gbif.py

# Rebuild (reads index.html + results/*.csv)
python scripts/build_combined_viewer.py
```

## Files

```
results/
  itp_fp_fits.csv        NIA ITP SR fits (64 curves)
  zims_all.csv           ZIMS zoo-animal SR fits (1 633 curves)
  zims_taxonomy.csv      GBIF taxonomy for all ZIMS species
scripts/
  fit_zims_one.py        Cluster fit script (one species × sex)
  run_zims_array.sh      LSF array-job wrapper
  build_combined_viewer.py  Builds viewer.html from index.html + results
  fetch_taxonomy_gbif.py    Queries GBIF for genus/family/order
  build_itp_overlay.py   ITP-only overlay builder (reference)
figures/
  zims_fig1_rms.png      RMS distribution by class
  zims_fig2_s_vs_L.png   Sharpness vs lifespan
  zims_fig3_identifiable_space.png  Identifiable parameter space
  zims_fig4_metadata.png Arm/ndims breakdown
```

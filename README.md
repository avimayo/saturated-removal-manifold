# Saturated-removal lifetime manifold — interactive viewer

An interactive 3-D view of the lifetime **constraint surface** of the saturated-removal aging model,
with fitted species placed in the dimensionless coordinates — removal $\rho_\beta$, production
$\rho_\eta$, diffusion $\rho_\epsilon$.

## ▶ Open the viewer

- **GitHub Pages:** https://avimayo.github.io/saturated-removal-manifold/
- **or via githack:** https://raw.githack.com/avimayo/saturated-removal-manifold/main/index.html

## What you're looking at

- The tan **surface** is the mean-lifetime constraint — a **two-dimensional manifold** in the three rates.
  Requiring the model's mean first-passage lifetime to equal the observed one gives one closed-form relation,
  $\tfrac12\rho_\eta-\rho_\beta = 1 - 1/\omega + e^{-\omega}/\omega$ with $\omega=(\tfrac12\rho_\eta-\rho_\beta)/\rho_\epsilon$.
  It depends on production and removal only through their net drift, so the rates are pinned to a surface
  (the production/removal split at fixed net drift stays free).
- **$\omega$** is the net-drift Péclet number — the lifetime-averaged net drift relative to diffusion — and marks
  position along the curve: $\omega\to+\infty$ deterministic (production-dominated), $\omega=0$ pure diffusion
  ($\rho_\epsilon=\tfrac12$), $\omega<0$ removal-dominated.
- The dark **ridge** is the **one-dimensional curve** the surface singles out in closed form,
  $(1-s-\rho_\epsilon)(\rho_\eta-\rho_\beta)=\tfrac12\rho_\eta\rho_\epsilon$ with $s=\tfrac12\rho_\eta-\rho_\beta$,
  labelled by $\omega$ (increasing toward the production end).
- The **points** are fitted species, colored by clade; they fall along the ridge.

Drag to rotate · scroll to zoom · use the legend to toggle layers, or the **show all / hide all data**
buttons to toggle every point at once · hover a point for its species and coordinates, or the ridge for $\omega$.

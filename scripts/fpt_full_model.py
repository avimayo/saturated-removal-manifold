"""First-passage statistics of the FULL saturated-removal model (NOT the high-damage reduction).

Dimensionless SDE (y=x/X_c, tau=t/L), with the genuine SATURATING, space-dependent removal:

    dy = a(y,tau) d tau + sqrt(2 rho_eps) dW,   a(y,tau) = rho_eta*tau - rho_beta * y/(y + kappa),
    kappa = k/X_c,   reflecting at y=0, absorbing at y=1.

Far from saturation (y <~ kappa) the removal becomes ~linear, rho_beta*y/kappa, vanishing at the wall
(an OU-like restoring force) -- qualitatively different from the reduced model's constant -rho_beta.
kappa -> 0 (k << X_c) recovers the reduced constant-removal model used elsewhere in the paper.

Method: forward Fokker-Planck, conservative Scharfetter-Gummel (exponential-fitting) flux so it is stable
at any local Peclet, integrated with BACKWARD EULER. (Crank-Nicolson is NOT L-stable: on the delta-function
start it oscillates, the p<0 clip then injects mass, and the survival comes out ~2x too large -> CV=0. BE
damps the high-frequency transient and conserves mass.) Survival S(tau)=int_0^1 p dy gives the moments by
mean=int S, <tau^2>=2 int tau S.

Validated: kappa->0 reproduces the exact constant-drift CV; the kappa-sweep matches a Langevin simulation of
the same full model to ~1%, at ~10x lower cost. Run `python3 scripts/fpt_full_model.py` for both checks.
"""
import numpy as np
from scipy.linalg import solve_banded
import mpmath as mp
mp.mp.dps = 30


def cv2_surrogate(om):
    """Exact constant-drift (reduced, rho_eta=0) CV^2 -- validation target for kappa->0."""
    om = mp.mpf(float(om))
    if abs(om) < 1e-9:
        return 2.0 / 3.0
    t1 = (om + mp.e**(-om) - 1) / om**2
    t2 = (om**2 - 4 + 2*mp.e**(-2*om) + 2*mp.e**(-om)*(3*om + 1)) / om**4
    return float(t2 / t1**2 - 1)


def _bernoulli(z):
    """B(z) = z/(e^z - 1), stable (series near 0, clipped exponent)."""
    z = np.clip(z, -500.0, 500.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(np.abs(z) < 1e-10, 1.0 - z/2.0, z/np.expm1(z))


def fpt_full(rho_eta, rho_beta, rho_eps, kappa, N=3000, dt=1e-3, tau_max=40.0, Gtol=1e-9,
             return_survival=False):
    """(mean, CV, S_end) of the first-passage time for the full saturating-removal model.

    If return_survival=True, returns (mean, CV, S_end, tau_grid, S_grid) instead, where
    S_grid is the survival G(tau)=int_0^1 p dy on the solver's own tau_grid (0..stop, step dt).
    The extra arrays are needed to fit the full empirical survival CURVE (see fp_curvefit_seed).
    """
    h = 1.0/N
    yf = (np.arange(N) + 0.5)*h          # cell-face midpoints y_{i+1/2}
    M = N; re = rho_eps; c = re/h**2

    def removal(yy):
        return rho_beta*yy/(yy + kappa)

    def build(a_face):                    # tridiagonal SG operator for drift sampled at the faces
        P = a_face*h/re
        Bp = _bernoulli(P); Bm = _bernoulli(-P)
        lo = np.zeros(M); di = np.zeros(M); up = np.zeros(M)
        i = np.arange(1, M-1)
        lo[i] = c*Bm[i-1]; di[i] = -c*(Bm[i] + Bp[i-1]); up[i] = c*Bp[i]
        di[0] = -c*Bm[0]; up[0] = c*Bp[0]                    # reflecting wall: zero flux at y=0
        lo[M-1] = c*Bm[M-2]; di[M-1] = -c*(Bm[M-1] + Bp[M-2])  # node N (y=1) absorbing, p_N=0
        return lo, di, up

    p = np.zeros(M); p[0] = 1.0/h         # delta at the reflecting wall
    nsteps = int(tau_max/dt)
    S = np.empty(nsteps+1); tl = np.empty(nsteps+1); S[0] = p.sum()*h; tl[0] = 0.0; last = 0
    for n in range(nsteps):
        lo1, di1, up1 = build(rho_eta*((n+1)*dt) - removal(yf))   # backward Euler: implicit at n+1
        ab = np.zeros((3, M)); ab[0, 1:] = -dt*up1[:-1]; ab[1, :] = 1.0 - dt*di1; ab[2, :-1] = -dt*lo1[1:]
        p = solve_banded((1, 1), ab, p); p[p < 0] = 0.0
        Sv = p.sum()*h; S[n+1] = Sv; tl[n+1] = (n+1)*dt; last = n+1
        if Sv < Gtol and (n+1)*dt > 0.3:
            break
    S = S[:last+1]; tl = tl[:last+1]
    mean = np.trapz(S, tl); m2 = 2*np.trapz(tl*S, tl)
    cv = np.sqrt(max(m2 - mean**2, 0.0))/mean
    if return_survival:
        return mean, cv, S[-1], tl, S
    return mean, cv, S[-1]


def fpt_full_survival(rho_eta, rho_beta, rho_eps, kappa, tau_eval, **kw):
    """Survival G(tau) of the full saturating-removal model on a requested tau-grid.

    Convenience wrapper around fpt_full(return_survival=True): solves once, then linearly
    interpolates the survival onto `tau_eval` (clamped to 1.0 at tau=0 and to the last solved
    value, ~0, beyond the absorbing tail). Returns (G_on_tau_eval, mean, cv, S_end).

    Pass solver controls (N, dt, tau_max, Gtol) through **kw. tau_max defaults high enough that
    the curve is fully absorbed; if a requested tau exceeds the integrated range the survival is
    held at its last (near-zero) value, which is the correct asymptote.
    """
    tau_eval = np.asarray(tau_eval, dtype=float)
    mean, cv, s_end, tg, Sg = fpt_full(rho_eta, rho_beta, rho_eps, kappa,
                                       return_survival=True, **kw)
    # left-clamp to S(0)=1, right-clamp to the last solved survival (the absorbed tail ~0)
    G = np.interp(tau_eval, tg, Sg, left=1.0, right=Sg[-1])
    return G, mean, cv, s_end


def mc_full(rho_eta, rho_beta, rho_eps, kappa, N=120000, tmax=18.0, seed=1):
    """Langevin Monte-Carlo of the same full model (cross-check)."""
    rng = np.random.default_rng(seed)
    dt = min(5e-4, 0.01/(2*rho_eps)); ns = int(tmax/dt)
    y = np.zeros(N); fpt = np.full(N, np.nan); act = np.ones(N, bool); sq = np.sqrt(2*rho_eps*dt)
    for n in range(ns):
        na = act.sum()
        if na == 0:
            break
        yy = y[act]
        y[act] = np.abs(yy + (rho_eta*(n*dt) - rho_beta*yy/(yy + kappa))*dt + sq*rng.standard_normal(na))
        cr = act & (y >= 1.0); fpt[cr] = (n+1)*dt; act &= ~cr
    f = fpt[~np.isnan(fpt)]
    return (f.mean(), f.std()/f.mean()) if len(f) > 100 else (np.nan, np.nan)


def main():
    print("== validation: kappa->0 recovers the reduced/exact constant-drift CV ==")
    m, cv, _ = fpt_full(0.0, 0.0, 0.5, 1e-6, tau_max=30.0)
    print(f"  pure diffusion : mean={m:.4f} (1.0)   CV={cv:.4f} ({np.sqrt(2/3):.4f})")
    for rb, re in [(1.0, 0.5), (2.0, 0.5), (4.0, 1.0)]:
        om = -rb/re; m, cv, _ = fpt_full(0.0, rb, re, 1e-6, tau_max=50.0)
        print(f"  omega={om:5.1f}    : mean={m:.4f}   CV_FP={cv:.4f}   CV_exact={np.sqrt(cv2_surrogate(om)):.4f}")
    print("\n== full model far from saturation: sweep kappa (rho_eta=3.85, rho_beta=2.47, rho_eps=0.435) ==")
    print(f"  {'kappa':>6s} {'mean_FP':>8s} {'CV_FP':>7s} | {'mean_MC':>8s} {'CV_MC':>7s}")
    for kappa in [1e-3, 0.05, 0.2, 0.5, 1.0]:
        mF, cF, _ = fpt_full(3.85, 2.47, 0.435, kappa, tau_max=18.0)
        mM, cM = mc_full(3.85, 2.47, 0.435, kappa)
        print(f"  {kappa:6.3f} {mF:8.4f} {cF:7.4f} | {mM:8.4f} {cM:7.4f}")


if __name__ == "__main__":
    main()

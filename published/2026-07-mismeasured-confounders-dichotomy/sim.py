"""
Simulation companion for
  "Adjusting for a mismeasured confounder is safe; adjusting for two is not:
   an exact residual-bias formula and a reliability-homogeneity dichotomy"
  (Chase's Journal, 2026-07)

What this checks / illustrates
------------------------------
Model (linear-Gaussian, scalar treatment X, k confounders U in R^k):
    U ~ N(0, Sigma_U)
    X = a'U + eps_X,          eps_X ~ N(0, s),            eps_X _|_ U
    Y = tau X + b'U + eps_Y,   eps_Y ~ N(0, v),            eps_Y _|_ (X,U)
    W = U + delta,             delta ~ N(0, Sigma_delta),  delta _|_ everything  (classical error)

Estimands (population OLS coefficient on X):
    naive : Y ~ X            bias  B_naive = (a' Sigma_U b) / (s + a' Sigma_U a)
    proxy : Y ~ (X, W)       bias  B_W     = (a' K b) / (s + a' K a),   K = Sigma_U Sigma_W^{-1} Sigma_delta
    oracle: Y ~ (X, U)       bias  0   (exact)
with Sigma_W = Sigma_U + Sigma_delta.

Claims verified numerically:
  (T1) closed form B_W matches the direct population block-inverse / OLS estimand,
       over random parameter draws (agreement ~1e-12), and matches finite-sample MC.
  (T2) SAFETY DICHOTOMY: |B_W| <= |B_naive| with matching sign for ALL (a,b,s)
       iff Sigma_delta ∝ Sigma_U (homogeneous reliability). We exhibit:
         - k=1 : always safe (F <= R), reproducing the scalar result;
         - k=2 homogeneous reliability : safe;
         - k=2 heterogeneous reliability : reversal, amplification, and
           MANUFACTURE (B_naive = 0 but B_W != 0).
  Figures: (1) manufactured-bias heatmap over (R1,R2); (2) overshoot-ratio curve
           showing the scalar ceiling 1 is breached (and sign flips) for k=2.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

SEED = 20260702
rng = np.random.default_rng(SEED)
FIGS = Path(__file__).resolve().parent / "figs"
FIGS.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# Population estimands (exact, via the true covariance matrices)
# ----------------------------------------------------------------------

def population_biases(Sigma_U, a, b, s, Sigma_delta):
    """Return (B_naive, B_W_closed, B_W_direct, B_oracle) as population estimands."""
    Sigma_U = np.asarray(Sigma_U, float)
    a = np.asarray(a, float); b = np.asarray(b, float)
    Sigma_delta = np.asarray(Sigma_delta, float)
    k = len(a)
    Sigma_W = Sigma_U + Sigma_delta

    A = a @ Sigma_U @ a + s
    B_naive = (a @ Sigma_U @ b) / A

    # K = Sigma_U Sigma_W^{-1} Sigma_delta
    K = Sigma_U @ np.linalg.solve(Sigma_W, Sigma_delta)
    B_W_closed = (a @ K @ b) / (s + a @ K @ a)

    # Direct: population OLS of Y on Z=(X,W). Bias = [M^{-1} g]_0, g = Cov(Z, b'U).
    cXW = Sigma_U @ a                       # Cov(X, W)
    M = np.block([[np.array([[A]]),      cXW[None, :]],
                  [cXW[:, None],         Sigma_W]])
    g = np.concatenate([[a @ Sigma_U @ b], Sigma_U @ b])
    B_W_direct = np.linalg.solve(M, g)[0]

    # Oracle: regress on (X,U); Sigma_delta=0 -> K=0 -> bias 0 exactly.
    B_oracle = 0.0
    return B_naive, B_W_closed, B_W_direct, B_oracle


def montecarlo_biases(Sigma_U, a, b, s, Sigma_delta, tau=2.0, n=200_000, reps=40, v=1.0):
    """Finite-sample OLS biases (mean over reps), for naive/proxy/oracle."""
    Sigma_U = np.asarray(Sigma_U, float); a = np.asarray(a, float); b = np.asarray(b, float)
    Sigma_delta = np.asarray(Sigma_delta, float); k = len(a)
    LU = np.linalg.cholesky(Sigma_U)
    Ld = np.linalg.cholesky(Sigma_delta + 1e-15 * np.eye(k))
    bn = bp = bo = 0.0
    for _ in range(reps):
        U = rng.standard_normal((n, k)) @ LU.T
        X = U @ a + np.sqrt(s) * rng.standard_normal(n)
        Y = tau * X + U @ b + np.sqrt(v) * rng.standard_normal(n)
        W = U + rng.standard_normal((n, k)) @ Ld.T
        # naive
        Zn = np.column_stack([np.ones(n), X])
        bn += np.linalg.lstsq(Zn, Y, rcond=None)[0][1]
        # proxy
        Zp = np.column_stack([np.ones(n), X, W])
        bp += np.linalg.lstsq(Zp, Y, rcond=None)[0][1]
        # oracle
        Zo = np.column_stack([np.ones(n), X, U])
        bo += np.linalg.lstsq(Zo, Y, rcond=None)[0][1]
    return (bn / reps - tau, bp / reps - tau, bo / reps - tau)


# ----------------------------------------------------------------------
# (T1) closed form == direct, over random draws
# ----------------------------------------------------------------------
def check_closed_form(ndraw=200_000):
    worst = 0.0
    for _ in range(ndraw):
        k = rng.integers(1, 5)
        # random SPD Sigma_U and diagonal Sigma_delta (classical error is per-coordinate)
        G = rng.standard_normal((k, k))
        Sigma_U = G @ G.T + 0.1 * np.eye(k)
        Sigma_delta = np.diag(rng.uniform(0.05, 3.0, k))
        a = rng.standard_normal(k)
        b = rng.standard_normal(k)
        s = rng.uniform(0.1, 3.0)
        _, cf, dr, _ = population_biases(Sigma_U, a, b, s, Sigma_delta)
        worst = max(worst, abs(cf - dr))
    return worst


# ----------------------------------------------------------------------
# (T2) safety dichotomy: search for violations under homogeneous vs
#      heterogeneous reliability
# ----------------------------------------------------------------------
def safety_scan(homogeneous, ndraw=300_000):
    """Return worst 'amplification' amount max(|B_W| - |B_naive|) and whether any
    sign reversal occurred, over random (a,b,s) for random covariance structure."""
    worst_amp = -np.inf
    reversal = False
    for _ in range(ndraw):
        k = 2
        G = rng.standard_normal((k, k))
        Sigma_U = G @ G.T + 0.1 * np.eye(k)
        if homogeneous:
            rho = rng.uniform(0.1, 3.0)
            Sigma_delta = rho * Sigma_U            # Sigma_delta ∝ Sigma_U
        else:
            Sigma_delta = np.diag(rng.uniform(0.05, 3.0, k))
        a = rng.standard_normal(k); b = rng.standard_normal(k)
        s = rng.uniform(0.1, 3.0)
        Bn, Bw, _, _ = population_biases(Sigma_U, a, b, s, Sigma_delta)
        worst_amp = max(worst_amp, abs(Bw) - abs(Bn))
        if Bn != 0 and np.sign(Bw) != np.sign(Bn) and abs(Bw) > 1e-9:
            reversal = True
    return worst_amp, reversal


# ----------------------------------------------------------------------
# Named scenarios for the results table
# ----------------------------------------------------------------------
def scenario_table():
    rows = []

    # (i) scalar baseline: one confounder, always safe (F <= R)
    rows.append(("scalar (k=1), R=0.5",
                 np.array([[1.0]]), np.array([1.0]), np.array([1.0]), 1.0,
                 np.diag([1.0])))  # Var U=1, w=1 -> R=0.5

    # (ii) k=2 homogeneous reliability: Sigma_delta ∝ Sigma_U -> safe
    Su = np.array([[1.0, 0.6], [0.6, 1.0]])
    rows.append(("k=2 homogeneous (Σδ=ΣU)",
                 Su, np.array([1.0, 0.8]), np.array([1.0, -0.5]), 1.0, 1.0 * Su))

    # (iii) k=2 MANUFACTURED bias: uncorrelated confounders, naive UNBIASED,
    #       heterogeneous reliability R1=0.9, R2=0.3 -> proxy biased
    #       Var U_j = 1 => Sigma_delta_jj = 1/R_j - 1
    R1, R2 = 0.9, 0.3
    Sd = np.diag([1 / R1 - 1, 1 / R2 - 1])
    rows.append(("k=2 manufacture (R=.9,.3)",
                 np.eye(2), np.array([1.0, 1.0]), np.array([1.0, -1.0]), 1.0, Sd))

    # (iv) k=2 SIGN REVERSAL: naive positive, proxy negative
    R1, R2 = 0.9, 0.2
    Sd = np.diag([1 / R1 - 1, 1 / R2 - 1])
    rows.append(("k=2 reversal (R=.9,.2)",
                 np.eye(2), np.array([1.0, 1.0]), np.array([1.0, -0.5]), 1.0, Sd))

    # (v) k=2 AMPLIFICATION (same sign, larger magnitude): confounder biases nearly
    #     cancel (a1b1=+2, a2b2=-1.9 -> small naive bias); the strongly biasing but
    #     poorly measured confounder 1 (R1=0.1) is barely attenuated, the offsetting
    #     well-measured confounder 2 (R2=0.9) is nearly removed -> |B_W| >> |B_naive|.
    R1, R2 = 0.1, 0.9
    Sd = np.diag([1 / R1 - 1, 1 / R2 - 1])
    rows.append(("k=2 amplify (R=.1,.9)",
                 np.eye(2), np.array([2.0, 1.9]), np.array([1.0, -1.0]), 1.0, Sd))

    out = []
    for name, Su, a, b, s, Sd in rows:
        Bn, Bw, Bwd, Bo = population_biases(Su, a, b, s, Sd)
        mBn, mBw, mBo = montecarlo_biases(Su, a, b, s, Sd)
        out.append(dict(name=name, Bn=Bn, Bw=Bw, Bwd=Bwd, Bo=Bo,
                        mBn=mBn, mBw=mBw, mBo=mBo))
    return out


# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def fig_manufacture_heatmap():
    """Two uncorrelated confounders with a1*b1 = +1, a2*b2 = -1 (naive bias = 0).
    Color = proxy-adjusted bias B_W over the reliability grid (R1,R2). Diagonal
    R1=R2 is the only safe line (B_W=0); off it, adjustment manufactures bias."""
    a = np.array([1.0, 1.0]); b = np.array([1.0, -1.0]); s = 1.0
    grid = np.linspace(0.05, 0.98, 200)
    Z = np.zeros((len(grid), len(grid)))
    for i, R2 in enumerate(grid):
        for j, R1 in enumerate(grid):
            Sd = np.diag([1 / R1 - 1, 1 / R2 - 1])
            Bn, Bw, _, _ = population_biases(np.eye(2), a, b, s, Sd)
            Z[i, j] = Bw
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    vmax = np.abs(Z).max()
    im = ax.pcolormesh(grid, grid, Z, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    ax.plot([0.05, 0.98], [0.05, 0.98], "k--", lw=1.4, label="$R_1=R_2$ (safe: $B_W=0$)")
    cs = ax.contour(grid, grid, Z, levels=[-0.3, -0.15, 0.15, 0.3], colors="k", linewidths=0.6)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.2f")
    ax.set_xlabel("reliability of proxy 1,  $R_1$")
    ax.set_ylabel("reliability of proxy 2,  $R_2$")
    ax.set_title("Bias manufactured by adjustment\n(naive is UNBIASED: $B_{\\rm naive}=0$ everywhere)")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    cb = fig.colorbar(im, ax=ax); cb.set_label("proxy-adjusted bias  $B_W$")
    fig.tight_layout()
    fig.savefig(FIGS / "manufacture_heatmap.png", dpi=150)
    plt.close(fig)


def fig_overshoot_ratio():
    """Fix confounder 1 well measured (R1=0.9); sweep R2 down. Plot B_W/B_naive.
    The scalar result confines this ratio to (0,1) (the shaded 'scalar-safe' band);
    with a 2nd confounder the curve leaves the band -- crossing 0 (reversal) and
    dropping below -1 (amplification: |B_W|>|B_naive|)."""
    a = np.array([1.0, 1.2]); b = np.array([1.0, -0.55]); s = 1.0  # naive bias > 0
    R1 = 0.9
    R2s = np.linspace(0.98, 0.05, 300)
    ratio = np.array([population_biases(np.eye(2), a, b, s,
                        np.diag([1/R1-1, 1/R2-1]))[1] /
                      population_biases(np.eye(2), a, b, s,
                        np.diag([1/R1-1, 1/R2-1]))[0] for R2 in R2s])
    fig, ax = plt.subplots(figsize=(5.8, 4.1))
    # scalar-safe band 0 < ratio < 1
    ax.axhspan(0.0, 1.0, color="#dfe7d8", alpha=0.8, zorder=0,
               label="scalar-safe band $0<B_W/B_{\\rm naive}<1$")
    ax.plot(R2s, ratio, lw=2.2, color="#b2182b", zorder=3)
    ax.axhline(1.0, color="gray", ls="--", lw=1.0)
    ax.axhline(0.0, color="gray", ls="--", lw=1.0)
    ax.axhline(-1.0, color="gray", ls=":", lw=1.0)
    # mark reversal onset (ratio crosses 0) and amplification onset (ratio < -1)
    i0 = int(np.argmin(np.abs(ratio)))
    ax.axvline(R2s[i0], color="#2166ac", ls=":", lw=1.2)
    ax.text(R2s[i0] + 0.005, -1.75, "sign reversal", color="#2166ac",
            rotation=90, va="center", fontsize=8)
    below = np.where(ratio < -1.0)[0]
    if len(below):
        ax.axvline(R2s[below[0]], color="#8c510a", ls=":", lw=1.2)
        ax.text(R2s[below[0]] + 0.005, -1.75, "amplification onset", color="#8c510a",
                rotation=90, va="center", fontsize=8)
    ax.set_xlabel("reliability of proxy 2,  $R_2$   (proxy 1 fixed at $R_1=0.9$)")
    ax.set_ylabel("$B_W\\,/\\,B_{\\rm naive}$")
    ax.set_title("Adjustment leaves the scalar-safe band with two confounders")
    ax.set_xlim(1.0, 0.05); ax.set_ylim(-2.2, 2.2)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(FIGS / "overshoot_ratio.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------
if __name__ == "__main__":
    print(f"seed = {SEED}\n")

    print("=== (T1) closed form vs direct population estimand ===")
    worst = check_closed_form(ndraw=200_000)
    print(f"  max |B_W(closed) - B_W(direct)| over 2e5 random draws (k=1..4): {worst:.2e}\n")

    print("=== (T2) safety dichotomy: worst amplification |B_W|-|B_naive| ===")
    amp_homo, rev_homo = safety_scan(homogeneous=True, ndraw=300_000)
    amp_het, rev_het = safety_scan(homogeneous=False, ndraw=300_000)
    print(f"  HOMOGENEOUS reliability (Σδ ∝ ΣU):  worst |B_W|-|B_naive| = {amp_homo:+.2e},  any sign reversal? {rev_homo}")
    print(f"  HETEROGENEOUS reliability        :  worst |B_W|-|B_naive| = {amp_het:+.3f},  any sign reversal? {rev_het}\n")

    print("=== scenario table (population theory vs Monte Carlo) ===")
    hdr = f"{'scenario':<28} {'B_naive':>16} {'B_W':>16} {'B_oracle':>14}"
    print(hdr); print("-" * len(hdr))
    for r in scenario_table():
        print(f"{r['name']:<28} "
              f"{r['Bn']:+8.4f}/{r['mBn']:+7.4f}  "
              f"{r['Bw']:+8.4f}/{r['mBw']:+7.4f}  "
              f"{r['Bo']:+6.3f}/{r['mBo']:+6.3f}")
    print("  (each cell: theory / Monte-Carlo)\n")

    print("=== figures ===")
    fig_manufacture_heatmap(); print("  wrote figs/manufacture_heatmap.png")
    fig_overshoot_ratio();     print("  wrote figs/overshoot_ratio.png")
    print("\ndone.")

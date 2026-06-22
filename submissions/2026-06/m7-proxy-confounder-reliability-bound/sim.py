"""
Adjusting for a noisy proxy of a confounder: how much bias do you actually remove?

Linear-Gaussian structural model
    U  ~ N(0, u)                      # unmeasured confounder
    X  = a*U + eps_X,  eps_X ~ N(0, s)   # treatment
    Y  = tau*X + b*U + eps_Y           # outcome (tau = causal effect of interest)
    W  = U + eps_W,    eps_W ~ N(0, w)   # classical-error proxy of U

We compare three population OLS estimands of the X-coefficient:
    naive   : Y ~ X            -> tau + B_full
    proxy   : Y ~ X + W        -> tau + B_W
    oracle  : Y ~ X + U        -> tau   (unbiased)

Theory (derived in note.md):
    A      = Var(X)            = a^2 u + s
    R      = reliability of W  = u / (u + w)
    Rx2    = a^2 u / A         = frac. of treatment variance explained by U
    B_full = b a u / A
    B_W    = b a u w / (a^2 u w + s (u + w))
    F = 1 - B_W/B_full = s u / (s u + w A)
      = R (1 - Rx2) / ( R (1 - Rx2) + (1 - R) )      <-- the fraction of bias removed
    Bound:  F <= R, with equality iff Rx2 = 0.

This script (1) checks the closed forms against Monte-Carlo OLS, and
(2) makes the two figures used in the note.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

rng = np.random.default_rng(20260622)
FIGS = Path(__file__).resolve().parent / "figs"
FIGS.mkdir(exist_ok=True)


# ----------------------------------------------------------------------
# Closed-form quantities
# ----------------------------------------------------------------------
def closed_form(u, s, a, b, w):
    A = a * a * u + s                       # Var(X)
    R = u / (u + w)                         # reliability of the proxy
    Rx2 = (a * a * u) / A                   # frac of Var(X) explained by U
    B_full = b * a * u / A
    B_W = b * a * u * w / (a * a * u * w + s * (u + w))
    F = 1.0 - (B_W / B_full if B_full != 0 else 0.0)
    F_alt = R * (1 - Rx2) / (R * (1 - Rx2) + (1 - R))
    return dict(A=A, R=R, Rx2=Rx2, B_full=B_full, B_W=B_W, F=F, F_alt=F_alt)


# ----------------------------------------------------------------------
# Monte-Carlo OLS estimands
# ----------------------------------------------------------------------
def ols_coef_on_x(Y, cols):
    """Population OLS: regress Y on [1, X, *cols]; return coefficient on X (index 1)."""
    n = len(Y)
    Xmat = np.column_stack([np.ones(n)] + cols)
    beta, *_ = np.linalg.lstsq(Xmat, Y, rcond=None)
    return beta[1]


def simulate(u, s, a, b, w, tau, v, n, reps):
    naive, proxy, oracle = [], [], []
    for _ in range(reps):
        U = rng.normal(0, np.sqrt(u), n)
        X = a * U + rng.normal(0, np.sqrt(s), n)
        Y = tau * X + b * U + rng.normal(0, np.sqrt(v), n)
        W = U + rng.normal(0, np.sqrt(w), n)
        naive.append(ols_coef_on_x(Y, [X]))
        proxy.append(ols_coef_on_x(Y, [X, W]))
        oracle.append(ols_coef_on_x(Y, [X, U]))
    return map(np.array, (naive, proxy, oracle))


# ----------------------------------------------------------------------
# (1) Numerical check of the closed forms
# ----------------------------------------------------------------------
print("=" * 72)
print("Monte-Carlo check of the closed-form bias (tau = 1.0 throughout)")
print("=" * 72)
tau, v, n, reps = 1.0, 1.0, 40000, 200

scenarios = [
    dict(u=1.0, s=1.0, a=1.0, b=1.0, w=0.5),   # moderate everything
    dict(u=1.0, s=0.25, a=2.0, b=1.5, w=1.0),  # U drives X strongly (high Rx2)
    dict(u=2.0, s=3.0, a=0.5, b=1.0, w=0.3),   # U weakly drives X (low Rx2)
]

for sc in scenarios:
    cf = closed_form(sc["u"], sc["s"], sc["a"], sc["b"], sc["w"])
    naive, proxy, oracle = simulate(tau=tau, v=v, n=n, reps=reps, **sc)
    mc_Bfull = naive.mean() - tau
    mc_BW = proxy.mean() - tau
    mc_oracle_bias = oracle.mean() - tau
    se = lambda arr: arr.std(ddof=1) / np.sqrt(reps)
    mc_F = 1 - mc_BW / mc_Bfull
    print(f"\nparams {sc}")
    print(f"  R = {cf['R']:.4f}   Rx2 = {cf['Rx2']:.4f}")
    print(f"  B_full : theory {cf['B_full']:+.4f}   MC {mc_Bfull:+.4f} (se {se(naive):.4f})")
    print(f"  B_W    : theory {cf['B_W']:+.4f}   MC {mc_BW:+.4f} (se {se(proxy):.4f})")
    print(f"  oracle : theory +0.0000   MC {mc_oracle_bias:+.4f} (se {se(oracle):.4f})")
    print(f"  F      : theory {cf['F']:.4f} (=alt {cf['F_alt']:.4f})   MC {mc_F:.4f}")
    print(f"  F <= R ? {cf['F']:.4f} <= {cf['R']:.4f}  -> {cf['F'] <= cf['R'] + 1e-12}")

# Cross-check F == F_alt over a random grid
print("\n" + "=" * 72)
print("Identity check  F == R(1-Rx2)/(R(1-Rx2)+1-R)  and bound F <= R over 100000 random params")
print("=" * 72)
M = 100_000
uu = rng.uniform(0.1, 5, M); ss = rng.uniform(0.1, 5, M)
aa = rng.uniform(-3, 3, M); bb = rng.uniform(-3, 3, M); ww = rng.uniform(0.01, 5, M)
A = aa**2 * uu + ss
R = uu / (uu + ww)
Rx2 = aa**2 * uu / A
F = ss * uu / (ss * uu + ww * A)
F_alt = R * (1 - Rx2) / (R * (1 - Rx2) + (1 - R))
print(f"  max |F - F_alt|            = {np.max(np.abs(F - F_alt)):.2e}")
print(f"  max (F - R)  (should be<=0)= {np.max(F - R):+.2e}")
print(f"  fraction with F <= R       = {np.mean(F <= R + 1e-12):.4f}")
print(f"  min F, max F               = {F.min():.4f}, {F.max():.4f}")


# ----------------------------------------------------------------------
# (2a) Figure: fraction removed F vs reliability R, several Rx2
# ----------------------------------------------------------------------
Rgrid = np.linspace(0.001, 0.999, 400)
fig, ax = plt.subplots(figsize=(6.2, 4.6))
ax.plot([0, 1], [0, 1], "k--", lw=1, label=r"naive expectation $F=R$")
colors = plt.cm.viridis(np.linspace(0.15, 0.85, 4))
for Rx2, c in zip([0.0, 0.3, 0.6, 0.9], colors):
    Fc = Rgrid * (1 - Rx2) / (Rgrid * (1 - Rx2) + (1 - Rgrid))
    ax.plot(Rgrid, Fc, color=c, lw=2, label=fr"$R_X^2={Rx2:.1f}$")

# overlay Monte-Carlo points (one per scenario) to show the curves are real
for sc in scenarios:
    cf = closed_form(sc["u"], sc["s"], sc["a"], sc["b"], sc["w"])
    naive, proxy, oracle = simulate(tau=tau, v=v, n=n, reps=reps, **sc)
    mc_F = 1 - (proxy.mean() - tau) / (naive.mean() - tau)
    ax.plot(cf["R"], mc_F, "o", ms=8, mfc="white", mec="crimson", mew=1.6, zorder=5)
ax.plot([], [], "o", ms=8, mfc="white", mec="crimson", mew=1.6, label="Monte-Carlo")

ax.set_xlabel("reliability of the proxy   $R = \\mathrm{Var}(U)/\\mathrm{Var}(W)$")
ax.set_ylabel("fraction of confounding bias removed   $F$")
ax.set_title("A proxy of reliability $R$ removes at most fraction $R$ of the bias")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.legend(loc="upper left", fontsize=9, framealpha=0.95)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGS / "fraction_removed.png", dpi=150)
print(f"\nwrote {FIGS/'fraction_removed.png'}")


# ----------------------------------------------------------------------
# (2b) Figure: the gap (R - F) as a heatmap over (R, Rx2)
# ----------------------------------------------------------------------
Rg = np.linspace(0.01, 0.99, 300)
Xg = np.linspace(0.0, 0.97, 300)
RR, XX = np.meshgrid(Rg, Xg)
FF = RR * (1 - XX) / (RR * (1 - XX) + (1 - RR))
GAP = RR - FF                                  # how much your proxy underperforms its reliability

fig, ax = plt.subplots(figsize=(6.2, 4.6))
pcm = ax.pcolormesh(RR, XX, GAP, shading="auto", cmap="magma_r")
cs = ax.contour(RR, XX, FF, levels=[0.1, 0.25, 0.5, 0.75, 0.9],
                colors="white", linewidths=1.0)
ax.clabel(cs, fmt=lambda x: f"F={x:.2f}", fontsize=8)
cb = fig.colorbar(pcm, ax=ax)
cb.set_label("shortfall  $R - F$  (bias you expected to remove but didn't)")
ax.set_xlabel("reliability of the proxy   $R$")
ax.set_ylabel("treatment variance explained by confounder   $R_X^2$")
ax.set_title("The stronger $U$ drives treatment, the worse the proxy performs")
fig.tight_layout()
fig.savefig(FIGS / "shortfall_heatmap.png", dpi=150)
print(f"wrote {FIGS/'shortfall_heatmap.png'}")
print("\nDONE.")

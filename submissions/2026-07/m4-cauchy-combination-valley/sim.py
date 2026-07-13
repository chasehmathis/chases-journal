"""
The robustness valley of the Cauchy combination test (ACAT).

Global-null simulation under a one-factor equicorrelated Gaussian copula.
We verify:
  (E)  exactness of the standard-Cauchy null at rho = 0 (independence) and
       rho = 1 (comonotonicity), for any finite m and equal weights;
  (V)  the "valley": the fixed-level size distortion S_m(rho, alpha) - alpha is
       zero at both endpoints and strictly positive (anti-conservative) in the
       interior, with an interior worst-case rho*;
  (S)  the exact symmetry T =d= -T for every rho (median stays 0);
  (A)  the tail-robustness of Liu & Xie (2020): the valley shrinks as alpha -> 0.

Outputs: figs/valley.png, figs/endpoints.png, and a printed summary that the
note quotes verbatim.
"""

import numpy as np
from scipy import stats

RNG = np.random.default_rng(20260713)


def cot(x):
    return 1.0 / np.tan(x)


def acat_stat(Z, w):
    """ACAT statistic from one-sided p-values p = 1 - Phi(Z).
    C_i = tan(pi(1/2 - p_i)) = cot(pi p_i); T = sum_i w_i C_i."""
    p = stats.norm.sf(Z)                     # 1 - Phi(Z), one-sided
    C = cot(np.pi * p)                        # exactly standard Cauchy under U(0,1)
    return C @ w


def simulate_T(rho, m, B, w):
    """B draws of the ACAT statistic under the global null, one-factor model:
       Z_i = sqrt(rho) F + sqrt(1-rho) eps_i,  F, eps_i ~ N(0,1) iid."""
    F = RNG.standard_normal((B, 1))
    eps = RNG.standard_normal((B, m))
    Z = np.sqrt(rho) * F + np.sqrt(max(0.0, 1.0 - rho)) * eps
    p = stats.norm.sf(Z)
    C = cot(np.pi * p)
    return C @ w


def size_at(rho, m, B, alpha_list, w):
    """Actual size P(T >= cot(pi alpha)) for each nominal alpha."""
    T = simulate_T(rho, m, B, w)
    return {a: float(np.mean(T >= cot(np.pi * a))) for a in alpha_list}, T


def mc_se(phat, B):
    return np.sqrt(phat * (1 - phat) / B)


# ---------------------------------------------------------------------------
# (E) + (S)  Endpoint exactness and symmetry
# ---------------------------------------------------------------------------
print("=" * 72)
print("(E) ENDPOINT EXACTNESS  (equal weights w_i = 1/m)")
print("=" * 72)
B_end = 4_000_000
alphas = [0.05, 0.01, 0.001]
for m in [2, 5, 20]:
    w = np.full(m, 1.0 / m)
    for rho, name in [(0.0, "independence rho=0"), (1.0, "comonotone   rho=1")]:
        s, T = size_at(rho, m, B_end, alphas, w)
        med = float(np.median(T))
        line = f"  m={m:3d}  {name}:  "
        line += "  ".join(
            f"a={a}: size={s[a]:.5f} (SE {mc_se(s[a], B_end):.5f})" for a in alphas
        )
        print(line)
    print(f"        median(T) at rho=0: {float(np.median(simulate_T(0.0, m, 500000, w))):+.4f}"
          f"   at rho=0.5: {float(np.median(simulate_T(0.5, m, 500000, w))):+.4f}"
          f"   at rho=1: {float(np.median(simulate_T(1.0, m, 500000, w))):+.4f}   (should be ~0)")

# KS check that T is standard Cauchy at the two endpoints (m=5)
w5 = np.full(5, 0.2)
for rho in [0.0, 1.0]:
    T = simulate_T(rho, 5, 200000, w5)
    ks = stats.kstest(T, 'cauchy')
    print(f"  KS test vs standard Cauchy, m=5, rho={rho}:  D={ks.statistic:.4f}, p={ks.pvalue:.3f}")

# ---------------------------------------------------------------------------
# (V) + (A)  The valley: size vs rho
# ---------------------------------------------------------------------------
print("=" * 72)
print("(V) THE VALLEY: actual size S_m(rho, alpha) across rho")
print("=" * 72)
rhos = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0])
B_val = 2_000_000
alpha_valley = [0.05, 0.01]
ms = [2, 5, 20, 100]

curves = {}      # (m, alpha) -> array of sizes over rhos
worst = {}       # (m, alpha) -> (rho*, size*)
for m in ms:
    w = np.full(m, 1.0 / m)
    for a in alpha_valley:
        curves[(m, a)] = np.zeros(len(rhos))
    for j, rho in enumerate(rhos):
        s, _ = size_at(rho, m, B_val, alpha_valley, w)
        for a in alpha_valley:
            curves[(m, a)][j] = s[a]
    for a in alpha_valley:
        c = curves[(m, a)]
        k = int(np.argmax(c))
        worst[(m, a)] = (rhos[k], c[k])
        print(f"  m={m:3d} alpha={a}:  size(0)={c[0]:.5f}  size(1)={c[-1]:.5f}  "
              f"max size={c[k]:.5f} at rho*={rhos[k]:.2f}  "
              f"(inflation x{c[k]/a:.2f})")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Figure 1: the valley (size vs rho), alpha = 0.05
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, a in zip(axes, alpha_valley):
    for m in ms:
        ax.plot(rhos, curves[(m, a)], marker="o", ms=3, label=f"m={m}")
    ax.axhline(a, color="k", ls="--", lw=1, label="nominal")
    ax.set_xlabel(r"latent correlation $\rho$")
    ax.set_ylabel("actual size")
    ax.set_title(f"Nominal $\\alpha={a}$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
fig.suptitle("The robustness valley: ACAT size is exact at $\\rho=0$ and $\\rho=1$, "
             "worst in between", fontsize=11)
fig.tight_layout()
fig.savefig("figs/valley.png", dpi=130)
plt.close(fig)

# Figure 2: endpoint exactness — empirical vs standard-Cauchy tail, and the
# interior departure, for m=20.
m = 20
w = np.full(m, 1.0 / m)
xs = np.linspace(0.5, 6, 200)                 # thresholds c
tail_cauchy = 0.5 - np.arctan(xs) / np.pi     # P(standard Cauchy > c)
fig, ax = plt.subplots(figsize=(6.2, 4.4))
ax.plot(xs, tail_cauchy, "k-", lw=2, label="standard Cauchy (ref.)")
for rho, style in [(0.0, "C0o"), (0.5, "C1s"), (0.9, "C3^"), (1.0, "C2d")]:
    T = simulate_T(rho, m, 1_000_000, w)
    emp = np.array([np.mean(T > c) for c in xs])
    ax.plot(xs, emp, style, ms=3, alpha=0.7, label=f"empirical, $\\rho={rho}$")
ax.set_yscale("log")
ax.set_xlabel("threshold $c$")
ax.set_ylabel(r"$P(T > c)$")
ax.set_title(f"Upper tail of the ACAT null ($m={m}$): exact at $\\rho\\in\\{{0,1\\}}$")
ax.legend(fontsize=8)
ax.grid(alpha=0.25, which="both")
fig.tight_layout()
fig.savefig("figs/endpoints.png", dpi=130)
plt.close(fig)

print("=" * 72)
print("Figures written: figs/valley.png, figs/endpoints.png")
print("=" * 72)

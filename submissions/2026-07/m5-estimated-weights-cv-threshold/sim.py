"""
When not to weight: the estimation penalty of plug-in inverse-variance weighting,
and a coefficient-of-variation threshold for when to just average.

Reproduces the two figures and all reported numbers in note.md.

Model: k independent groups; group i gives n i.i.d. Gaussian observations with
unknown variance sigma_i^2 and common mean theta (WLOG theta = 0). We compare
three combined estimators of theta built from the group means xbar_i:

  - EQUAL   : (1/k) sum_i xbar_i
  - PLUGIN  : sum_i what_i xbar_i,  what_i propto 1/s_i^2   (Graybill-Deal / feasible IVW)
  - ORACLE  : sum_i w_i^* xbar_i,   w_i^* propto 1/sigma_i^2 (known variances)

Part A (homoscedastic, sigma_i^2 == sigma^2): verifies
    Var(PLUGIN)/Var(EQUAL) = k * E[sum what_i^2] = 1 + 2(k-1)/(k(n-1)) + O(1/n^2).

Part B (heteroscedastic): verifies the crossover
    PLUGIN beats EQUAL  <=>  c_v^2 := Var(sigma_i^2)/mean(sigma_i^2)^2  >  2(k-1)/(k(n-1)).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

SEED = 20260716
rng = np.random.default_rng(SEED)
FIGS = "figs"

import os
os.makedirs(FIGS, exist_ok=True)


# ----------------------------------------------------------------------------
# Part A: the homoscedastic estimation penalty  k * E[sum what_i^2]
# ----------------------------------------------------------------------------
# With sigma_i^2 == sigma^2, the plug-in weights depend only on the sample
# variances s_i^2, which for Gaussian data are independent of the group means.
# Hence Var(PLUGIN) = sigma^2 * E[sum what_i^2], and the *ratio* to Var(EQUAL)
# = sigma^2/k is exactly  k * E[sum what_i^2], a pure number in (1, infinity).
# what_i = (1/s_i^2)/sum_j(1/s_j^2); the sigma^2 cancels, so we simulate
# s_i^2 proportional to chi^2_{n-1} directly.

def penalty_factor(k, n, reps):
    nu = n - 1
    V = rng.chisquare(nu, size=(reps, k))          # s_i^2 ~ (sigma^2/nu) chi^2_nu
    W = 1.0 / V
    w = W / W.sum(axis=1, keepdims=True)
    sw2 = (w ** 2).sum(axis=1)                      # sum_i what_i^2
    kR = k * sw2.mean()
    se = k * sw2.std() / np.sqrt(reps)              # Monte Carlo s.e. of the factor
    return kR, se


print("=" * 74)
print("PART A  Homoscedastic penalty  Var(PLUGIN)/Var(EQUAL) = k*E[sum what^2]")
print("        prediction: 1 + 2(k-1)/(k(n-1))")
print("=" * 74)
REPS_A = 1_500_000
tableA = []
for k in (2, 5, 10):
    for n in (6, 11, 21, 51, 101):
        kR, se = penalty_factor(k, n, REPS_A)
        pred = 1 + 2 * (k - 1) / (k * (n - 1))
        tableA.append((k, n, kR, se, pred))
        print(f"  k={k:2d} n={n:4d} | MC ratio = {kR:.4f} (+/-{se:.4f}) | "
              f"1+2(k-1)/(k(n-1)) = {pred:.4f}")

# Convergence of the leading coefficient: (n-1)*(ratio-1) -> 2(k-1)/k.
print("\n  Leading-coefficient check  (n-1)*(ratio-1) -> 2(k-1)/k:")
for k in (2, 5, 10):
    row = []
    for n in (51, 101, 201, 401):
        kR, _ = penalty_factor(k, n, REPS_A)
        row.append((n, (n - 1) * (kR - 1)))
    tail = "  ".join(f"n={n}:{val:.3f}" for n, val in row)
    print(f"    k={k:2d}: target 2(k-1)/k={2*(k-1)/k:.3f} | {tail}")


# --- Figure 1: penalty vs n, overlaid with the closed-form leading term --------
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ns = np.array([4, 6, 8, 11, 16, 21, 31, 51, 101])
colors = {2: "#1b6ca8", 5: "#c0392b", 10: "#2e7d32"}
for k in (2, 5, 10):
    mc = []
    for n in ns:
        kR, _ = penalty_factor(k, n, 1_000_000)
        mc.append(100 * (kR - 1))
    ax.plot(ns, mc, "o", color=colors[k], ms=6, label=f"$k={k}$ (simulated)")
    pred = 100 * 2 * (k - 1) / (k * (ns - 1))
    ax.plot(ns, pred, "-", color=colors[k], lw=1.6, alpha=0.8)
ax.set_xlabel("group size $n$")
ax.set_ylabel("excess variance of plug-in IVW over\nsimple averaging (%)")
ax.set_title("Homoscedastic penalty: markers = simulation, lines = $2(k{-}1)/(k(n{-}1))$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.grid(True, which="both", alpha=0.25)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{FIGS}/penalty.png", dpi=140); plt.close(fig)


# ----------------------------------------------------------------------------
# Part B: heteroscedastic crossover in the coefficient of variation
# ----------------------------------------------------------------------------
def variances_with_cv2(k, cv2):
    """Deterministic variance profile sigma_i^2 with mean 1 and squared
    coefficient of variation exactly cv2 (log-linear spread)."""
    if cv2 < 1e-9:
        return np.ones(k)
    z = np.linspace(-1, 1, k)
    f = lambda g: (np.exp(g * z).var() / np.exp(g * z).mean() ** 2) - cv2
    g = brentq(f, 1e-6, 40)
    v = np.exp(g * z)
    return v / v.mean()


def empirical_variances(k, n, sig2, reps):
    sig = np.sqrt(sig2)
    X = rng.standard_normal((reps, k, n)) * sig[None, :, None]
    xbar = X.mean(axis=2)
    s2 = X.var(axis=2, ddof=1)
    eq = xbar.mean(axis=1)
    wp = 1.0 / s2; wp /= wp.sum(axis=1, keepdims=True)
    plug = (wp * xbar).sum(axis=1)
    wo = (1.0 / sig2)[None, :]; wo = wo / wo.sum()
    orac = (wo * xbar).sum(axis=1)
    return eq.var(), plug.var(), orac.var()


print("\n" + "=" * 74)
print("PART B  Heteroscedastic crossover (equal vs plug-in IVW)")
print("=" * 74)
k, n = 8, 11
thresh = 2 * (k - 1) / (k * (n - 1))
print(f"  k={k}, n={n}: predicted crossover c_v^2 = 2(k-1)/(k(n-1)) = {thresh:.4f}")
print(f"  {'c_v^2':>7} {'H=eq/orc':>9} {'Var eq':>10} {'Var plug':>10} "
      f"{'Var orac':>10} {'plug<eq?':>9}")

REPS_B = 500_000
cv2_grid = [0.00, 0.02, 0.05, 0.08, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.35, 0.50]
rows = []
for cv2 in cv2_grid:
    sig2 = variances_with_cv2(k, cv2)
    H = sig2.mean() * np.mean(1.0 / sig2)
    ve, vp, vo = empirical_variances(k, n, sig2, REPS_B)
    rows.append((cv2, H, ve, vp, vo))
    print(f"  {cv2:7.3f} {H:9.4f} {ve:10.6f} {vp:10.6f} {vo:10.6f} "
          f"{str(vp < ve):>9}")

# locate empirical crossover (linear interpolation of vp-ve through 0)
d = np.array([vp - ve for (_, _, ve, vp, vo) in rows])
c = np.array([r[0] for r in rows])
xover = None
for i in range(len(d) - 1):
    if d[i] > 0 >= d[i + 1]:
        xover = c[i] + (c[i + 1] - c[i]) * d[i] / (d[i] - d[i + 1])
        break
print(f"\n  Empirical crossover c_v^2 ~= {xover:.3f}   (predicted {thresh:.3f})")


# --- Figure 2: crossover picture ----------------------------------------------
cv2_dense = np.linspace(0.0, 0.5, 26)
ve_l, vp_l, vo_l = [], [], []
for cv2 in cv2_dense:
    sig2 = variances_with_cv2(k, cv2)
    ve, vp, vo = empirical_variances(k, n, sig2, 250_000)
    ve_l.append(ve); vp_l.append(vp); vo_l.append(vo)
ve_l, vp_l, vo_l = map(np.array, (ve_l, vp_l, vo_l))

fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(cv2_dense, ve_l, "-", color="#555", lw=2, label="simple average (EQUAL)")
ax.plot(cv2_dense, vp_l, "-", color="#c0392b", lw=2, label="plug-in IVW (PLUGIN)")
ax.plot(cv2_dense, vo_l, "--", color="#1b6ca8", lw=1.6, label="oracle IVW (known var.)")
ax.axvline(thresh, color="k", ls=":", lw=1.4)
ax.annotate(f"predicted crossover\n$c_v^2 = 2(k{{-}}1)/(k(n{{-}}1)) = {thresh:.3f}$",
            xy=(thresh, ve_l[0]), xytext=(thresh + 0.03, ve_l[0] * 1.02),
            fontsize=9, va="center")
ax.set_xlabel(r"heterogeneity of true variances  $c_v^2 = \mathrm{Var}(\sigma_i^2)/\overline{\sigma^2}^2$")
ax.set_ylabel("variance of combined estimator")
ax.set_title(f"When to weight ($k={k}$, $n={n}$): plug-in IVW beats averaging only past the threshold")
ax.grid(True, alpha=0.25)
ax.legend(frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig(f"{FIGS}/crossover.png", dpi=140); plt.close(fig)

print("\nWrote figs/penalty.png and figs/crossover.png")
print(f"SEED = {SEED}")

"""
Sign vs. size in the significance filter.

Model: a two-sided z-test at level alpha, critical value c = Phi^{-1}(1 - alpha/2).
The observed standardized estimate is Z ~ N(d, 1) with true standardized effect
d = theta/SE > 0 (WLOG positive). Conditioning on "significant" means |Z| > c.

Two selection-induced error summaries (Gelman & Carlin 2014):
  Type-M (exaggeration ratio):  M(d) = E[|Z| | |Z|>c] / d
  Type-S (sign-error rate):     S(d) = P(Z < -c | |Z|>c)

This script:
  1. Verifies the closed forms for M(d) and S(d) against Monte Carlo.
  2. Verifies the low-power asymptotic law  M(d) ~ r(c)/d,  r(c)=phi(c)/Phi(-c).
  3. Locates the "reliability window": the band of power over which S(d) < 0.01
     while M(d) > 1.3, and reports its endpoints.
  4. Produces two figures.

Everything is deterministic given the seed.
"""

import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260727)

alpha = 0.05
c = norm.ppf(1 - alpha / 2)            # 1.959963...
phi = norm.pdf
Phi = norm.cdf


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------
def power(d):
    """Two-sided power at standardized effect d."""
    return Phi(d - c) + Phi(-d - c)

def typeM(d):
    """Exaggeration ratio M(d) = E[|Z| | |Z|>c] / d, closed form."""
    pw = power(d)
    num = d * (Phi(d - c) - Phi(-c - d)) + (phi(c - d) + phi(c + d))
    return (num / pw) / d

def typeS(d):
    """Sign-error rate S(d) = P(Z<-c | |Z|>c), closed form."""
    return Phi(-c - d) / power(d)

def cond_absZ_mean(d):
    """E[|Z| | |Z|>c], closed form (the numerator of Type-M times d)."""
    pw = power(d)
    return (d * (Phi(d - c) - Phi(-c - d)) + (phi(c - d) + phi(c + d))) / pw

def r_of_c(c):
    """Inverse Mills ratio at the critical value: phi(c)/Phi(-c)."""
    return phi(c) / Phi(-c)


# ---------------------------------------------------------------------------
# 1. Monte-Carlo check of the closed forms
# ---------------------------------------------------------------------------
print("=" * 70)
print(f"alpha = {alpha}, two-sided critical value c = {c:.6f}")
print(f"inverse Mills ratio at c:  r(c) = phi(c)/Phi(-c) = {r_of_c(c):.6f}")
print("=" * 70)

N = 20_000_000
print("\n[1] Monte-Carlo verification of closed forms (N = %d draws/effect):" % N)
print(f"{'d':>6} {'power':>8} {'M_exact':>9} {'M_MC':>9} {'S_exact':>10} {'S_MC':>11}")
for d in [0.5, 1.0, 1.5, 2.0, 2.8]:
    Z = rng.normal(d, 1.0, size=N)
    sig = np.abs(Z) > c
    M_mc = np.abs(Z[sig]).mean() / d
    S_mc = np.mean(Z[sig] < -c)
    print(f"{d:6.2f} {power(d):8.4f} {typeM(d):9.4f} {M_mc:9.4f} "
          f"{typeS(d):10.5f} {S_mc:11.6f}")


# ---------------------------------------------------------------------------
# 2. Low-power asymptotic  M(d) ~ r(c)/d
# ---------------------------------------------------------------------------
print("\n[2] Low-power asymptotic  M(d) * d  ->  r(c) = %.6f  as d -> 0:" % r_of_c(c))
print(f"{'d':>8} {'power':>9} {'M(d)*d':>10} {'ratio to r(c)':>14}")
for d in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
    Md = typeM(d) * d
    print(f"{d:8.3f} {power(d):9.4f} {Md:10.5f} {Md / r_of_c(c):14.5f}")


# ---------------------------------------------------------------------------
# 3. The reliability window:  S(d) < eps_S  AND  M(d) > 1 + eps_M
# ---------------------------------------------------------------------------
eps_S = 0.01     # sign-error rate below 1%
thr_M = 1.30     # magnitude inflated by more than 30%

# S is increasing as d decreases; find d where S(d) = eps_S.
d_lo = brentq(lambda d: typeS(d) - eps_S, 0.2, 6.0)       # S(d_lo) = eps_S
# M is decreasing in d; find d where M(d) = thr_M.
d_hi = brentq(lambda d: typeM(d) - thr_M, 0.2, 6.0)       # M(d_hi) = thr_M

print("\n[3] Reliability window  {S(d) < %.2f  and  M(d) > %.2f}:" % (eps_S, thr_M))
print(f"    lower edge (S = {eps_S}):  d = {d_lo:.4f},  power = {power(d_lo):.4f},"
      f"  M = {typeM(d_lo):.4f}")
print(f"    upper edge (M = {thr_M}):  d = {d_hi:.4f},  power = {power(d_hi):.4f},"
      f"  S = {typeS(d_hi):.6f}")
print(f"    => window in POWER: ({power(d_lo):.3f}, {power(d_hi):.3f})")

# A table of headline (power -> M, S) pairs.
print("\n    Headline pairs:")
print(f"{'power':>8} {'d':>7} {'M (exaggeration)':>18} {'S (sign error)':>16}")
for pw_target in [0.05001, 0.10, 0.20, 0.50, 0.80]:
    if pw_target <= alpha:
        d_star = 1e-6
    else:
        d_star = brentq(lambda d: power(d) - pw_target, 1e-6, 10.0)
    print(f"{power(d_star):8.3f} {d_star:7.3f} {typeM(d_star):18.3f} "
          f"{typeS(d_star):16.5f}")


# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------
dd = np.linspace(0.02, 4.0, 800)
pw = power(dd)
M = typeM(dd)
S = typeS(dd)

# Fig 1: M and S vs power, window shaded.
fig, ax1 = plt.subplots(figsize=(7.2, 4.4))
order = np.argsort(pw)
ax1.plot(pw[order], M[order], color="C3", lw=2.2, label="Type-M  (exaggeration)")
ax1.axhline(1.0, color="C3", ls=":", lw=1, alpha=0.6)
ax1.set_xlabel("power  $\\pi(d)$")
ax1.set_ylabel("exaggeration ratio  $M(d)$", color="C3")
ax1.tick_params(axis="y", labelcolor="C3")
ax1.set_ylim(0.8, 6.0)
ax1.set_xlim(alpha, 1.0)

ax2 = ax1.twinx()
ax2.plot(pw[order], S[order], color="C0", lw=2.2, label="Type-S  (sign error)")
ax2.set_ylabel("sign-error rate  $S(d)$  (log)", color="C0")
ax2.tick_params(axis="y", labelcolor="C0")
ax2.set_yscale("log")
ax2.set_ylim(1e-6, 1.0)

# Shade the reliability window (in power).
ax1.axvspan(power(d_lo), power(d_hi), color="grey", alpha=0.13, zorder=0)
ax1.text(0.5 * (power(d_lo) + power(d_hi)), 5.3,
         "reliability window\n(sign trusted, size inflated)",
         ha="center", va="top", fontsize=8.5, color="0.25")

ax1.set_title("Sign vs. size under the significance filter  ($\\alpha=0.05$, two-sided)")
lines = ax1.get_lines()[:1] + ax2.get_lines()[:1]
ax1.legend(lines, [l.get_label() for l in lines], loc="upper right", fontsize=9)
fig.tight_layout()
fig.savefig("figs/sign_vs_size.png", dpi=150)
plt.close(fig)

# Fig 2: low-power asymptotic  M(d)*d -> r(c).
fig, ax = plt.subplots(figsize=(7.2, 4.0))
dsmall = np.linspace(0.01, 2.0, 500)
ax.plot(dsmall, typeM(dsmall) * dsmall, color="C2", lw=2.2,
        label="$M(d)\\,d$  (exact)")
ax.axhline(r_of_c(c), color="k", ls="--", lw=1.4,
           label="$r(c)=\\varphi(c)/\\Phi(-c)=%.3f$" % r_of_c(c))
ax.set_xlabel("true standardized effect  $d$")
ax.set_ylabel("$M(d)\\,d$")
ax.set_title("Low-power law: the exaggeration blows up as $r(c)/d$")
ax.set_xlim(0, 2.0)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig("figs/asymptotic.png", dpi=150)
plt.close(fig)

print("\nWrote figs/sign_vs_size.png and figs/asymptotic.png")
print("Done.")

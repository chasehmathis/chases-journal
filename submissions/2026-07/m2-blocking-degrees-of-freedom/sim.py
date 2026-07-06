"""
Simulation companion for:
  "The degrees-of-freedom tax on blocking: a sharp threshold, a (z^2+1)/4n
   rule, and why finer blocking is the hardest to justify"
  Chase's Journal, 2026-07.

We study a completely randomized experiment on N = g*b units grouped into g
blocks of size b, with a one-way random-effects structure of intra-class
correlation (ICC) rho.  Two designs are compared on the SAME N units:

  * BLOCKED  : assign b/2 treated / b/2 control WITHIN each block; analyze with
               block fixed effects.  True Var(tau_hat) = 4 sigma^2 (1-rho) / N,
               residual df = N - g - 1.
  * CR (unblocked): assign N/2 treated / N/2 control completely at random,
               ignoring blocks; two-sample analysis.
               True Var(tau_hat) = 4 sigma^2 / N, residual df = N - 2.

Blocking cuts the variance by the factor (1-rho) but spends g-1 extra error df,
inflating the t critical value.  The break-even ICC for equal (oracle) CI
half-width is

        rho*(g) = 1 - ( t_{N-2, a/2} / t_{N-g-1, a/2} )^2 ,           (Prop 1)

with, for PAIRING (b=2, g=n, N=2n),
        rho*(n) = 1 - ( t_{2n-2} / t_{n-1} )^2  ~  (z^2 + 1) / (4(n-1)).  (Cor)

This script:
  (A) tabulates rho*(n) for pairing, exact vs the (z^2+1)/(4(n-1)) rule;
  (B) simulates true power & mean CI half-width of paired vs CR designs across
      rho, locates the empirical power/width crossings, and checks them against
      rho*;
  (C) shows rho*(g) is DECREASING in block size b (increasing in g): pairing is
      the most demanding blocking, i.e. the "block as finely as possible"
      instinct carries the largest df tax.  Also contrasts with the classical
      Fisher / Cochran-Cox variance-information relative-efficiency df factor.

Everything printed here is reproduced verbatim in note.md.
"""

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

rng = np.random.default_rng(20260706)
ALPHA = 0.05
Z = stats.norm.ppf(1 - ALPHA / 2)              # 1.959963...
FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
os.makedirs(FIGDIR, exist_ok=True)


def tcrit(df):
    return stats.t.ppf(1 - ALPHA / 2, df)


# ----------------------------------------------------------------------------
# (A) Pairing threshold: exact rho*(n) vs the (z^2+1)/(4(n-1)) rule of thumb
# ----------------------------------------------------------------------------
def rho_star_pairing_exact(n):
    return 1.0 - (tcrit(2 * n - 2) / tcrit(n - 1)) ** 2


def rho_star_pairing_rule(n):
    return (Z ** 2 + 1) / (4 * (n - 1))


print("=" * 72)
print("(A) Pairing break-even ICC  rho*(n)   [alpha = 0.05, z^2+1 = %.4f]" % (Z ** 2 + 1))
print("=" * 72)
print(f"{'n pairs':>8} {'exact rho*':>12} {'(z^2+1)/4(n-1)':>16} {'abs err':>10}")
ns_tab = [3, 4, 5, 8, 10, 15, 20, 30, 50, 100]
for n in ns_tab:
    ex = rho_star_pairing_exact(n)
    ru = rho_star_pairing_rule(n)
    print(f"{n:>8} {ex:>12.4f} {ru:>16.4f} {abs(ex - ru):>10.4f}")


# ----------------------------------------------------------------------------
# (B) Monte-Carlo: true power and mean CI half-width, paired vs CR, over rho.
#     Same N = 2n units both designs; fixed treatment effect tau.
# ----------------------------------------------------------------------------
def gen_units(n, rho, tau_on_treated, assign, reps):
    """
    Generate `reps` datasets of N = 2n units, ICC rho (unit var = 1):
      unit value = sqrt(rho)*block_effect + sqrt(1-rho)*noise,
      blocks of size 2 (n blocks).  `assign` in {'paired','cr'} sets W.
    Returns tau_hat, se_hat, df for each rep (vectorized where easy).
    """
    b_eff = rng.standard_normal((reps, n)) * np.sqrt(rho)          # per-pair
    noise = rng.standard_normal((reps, n, 2)) * np.sqrt(1 - rho)
    Y = b_eff[:, :, None] + noise                                  # (reps,n,2)

    if assign == "paired":
        # unit 0 treated, unit 1 control within each pair
        Y = Y.copy()
        Y[:, :, 0] += tau_on_treated
        D = Y[:, :, 0] - Y[:, :, 1]                                # (reps,n)
        tau_hat = D.mean(1)
        se_hat = D.std(1, ddof=1) / np.sqrt(n)
        df = n - 1
        return tau_hat, se_hat, df
    else:  # complete randomization over the 2n units, ignore pairs
        flat = Y.reshape(reps, 2 * n)                              # (reps,2n)
        # random half treated
        order = rng.random((reps, 2 * n)).argsort(1)
        treated = order[:, :n]
        control = order[:, n:]
        r = np.arange(reps)[:, None]
        Yt = flat[r, treated] + tau_on_treated
        Yc = flat[r, control]
        tau_hat = Yt.mean(1) - Yc.mean(1)
        vt = Yt.var(1, ddof=1)
        vc = Yc.var(1, ddof=1)
        sp2 = ((n - 1) * vt + (n - 1) * vc) / (2 * n - 2)          # pooled
        se_hat = np.sqrt(sp2 * (2.0 / n))
        df = 2 * n - 2
        return tau_hat, se_hat, df


def power_and_width(n, rho, tau, assign, reps):
    tau_hat, se_hat, df = gen_units(n, rho, tau, assign, reps)
    tcr = tcrit(df)
    tstat = tau_hat / se_hat
    reject = np.abs(tstat) > tcr
    half = tcr * se_hat
    return reject.mean(), half.mean()


REPS = 40000
n_demo = 8                     # 8 pairs => N = 16
# choose tau so CR power is ~0.6 at rho=0 (moderate)
true_var_cr = 4.0 / (2 * n_demo)          # = 2/n_demo (sigma^2=1)
tau = 0.62 * np.sqrt(true_var_cr) * (Z + stats.norm.ppf(0.60))   # rough targeting
rhos = np.linspace(0.0, 0.6, 25)

print("\n" + "=" * 72)
print(f"(B) Monte-Carlo power / half-width, n={n_demo} pairs (N={2*n_demo}), "
      f"tau={tau:.3f}, reps={REPS}")
print("=" * 72)
rho_star_demo = rho_star_pairing_exact(n_demo)
print(f"predicted rho*({n_demo}) = {rho_star_demo:.4f}  "
      f"[(z^2+1)/4(n-1) rule: {rho_star_pairing_rule(n_demo):.4f}]")

pw_p, pw_c, hw_p, hw_c = [], [], [], []
for r in rhos:
    pp, hp = power_and_width(n_demo, r, tau, "paired", REPS)
    pc, hc = power_and_width(n_demo, r, tau, "cr", REPS)
    pw_p.append(pp); pw_c.append(pc); hw_p.append(hp); hw_c.append(hc)
pw_p = np.array(pw_p); pw_c = np.array(pw_c)
hw_p = np.array(hw_p); hw_c = np.array(hw_c)


def crossing(x, ya, yb):
    """first x where ya-yb changes sign (linear interp)."""
    d = ya - yb
    for i in range(len(d) - 1):
        if d[i] == 0:
            return x[i]
        if d[i] * d[i + 1] < 0:
            t = d[i] / (d[i] - d[i + 1])
            return x[i] + t * (x[i + 1] - x[i])
    return np.nan


# power: paired beats CR when pw_p > pw_c  => crossing of (pw_p - pw_c)
x_pow = crossing(rhos, pw_p, pw_c)
# width: paired shorter when hw_p < hw_c  => crossing of (hw_c - hw_p)
x_wid = crossing(rhos, hw_c, hw_p)
print(f"empirical power  crossing (paired overtakes CR): rho = {x_pow:.4f}")
print(f"empirical width  crossing (paired shorter)     : rho = {x_wid:.4f}")
print(f"oracle threshold rho*                          : rho = {rho_star_demo:.4f}")

# Also print a couple of concrete rows near/around rho*
for target in [0.0, rho_star_demo, 0.4]:
    j = int(np.argmin(np.abs(rhos - target)))
    print(f"  rho={rhos[j]:.3f}:  power paired={pw_p[j]:.3f} CR={pw_c[j]:.3f} | "
          f"half-width paired={hw_p[j]:.3f} CR={hw_c[j]:.3f}")


# ----------------------------------------------------------------------------
# (C) Non-monotonicity in block size, and the Fisher/Cochran-Cox contrast.
#     Fix N; vary block size b (=> g = N/b blocks).  rho*(g) below.
# ----------------------------------------------------------------------------
def rho_star_block(N, g):                       # CI-width / power threshold
    return 1.0 - (tcrit(N - 2) / tcrit(N - g - 1)) ** 2


def rho_star_fisher(N, g):
    """Break-even from the classical Fisher / Cochran-Cox relative-efficiency
    df correction.  RE of blocked-to-CR = factor * (sigma^2_CR / sigma^2_blk)
    = factor / (1-rho), with the df factor

        factor = (f_b+1)(f_c+3) / [(f_b+3)(f_c+1)] < 1   (penalizes fewer df),

    f_b = N-g-1 (blocked error df), f_c = N-2 (CR error df).  RE = 1 gives
    1-rho = factor, i.e. rho* = 1 - factor."""
    fb, fc = N - g - 1, N - 2
    factor = (fb + 1) * (fc + 3) / ((fb + 3) * (fc + 1))
    return 1.0 - factor


N_fixed = 60
block_sizes = [2, 3, 4, 5, 6, 10, 12, 15, 20, 30]
print("\n" + "=" * 72)
print(f"(C) Fixed N={N_fixed}: break-even ICC vs block size b  (g=N/b blocks)")
print("=" * 72)
print(f"{'b':>4} {'g':>5} {'err df':>7} {'rho* (width/power)':>19} {'rho* (Fisher RE)':>17}")
bs_ok, rs_width, rs_fisher = [], [], []
for b in block_sizes:
    if N_fixed % b:
        continue
    g = N_fixed // b
    rw = rho_star_block(N_fixed, g)
    rf = rho_star_fisher(N_fixed, g)
    bs_ok.append(b); rs_width.append(rw); rs_fisher.append(rf)
    print(f"{b:>4} {g:>5} {N_fixed-g-1:>7} {rw:>19.4f} {rf:>17.4f}")
print("=> rho*(width/power) is DECREASING in b: pairing (b=2) is the hardest to justify.")
print("=> Fisher/Cochran-Cox RE threshold is uniformly LOWER (less demanding).")

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
# Fig 1: rho*(n) pairing exact vs rule
ns = np.arange(3, 101)
ex = np.array([rho_star_pairing_exact(n) for n in ns])
ru = np.array([rho_star_pairing_rule(n) for n in ns])
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.plot(ns, ex, "-", lw=2, label=r"exact  $1-(t_{2n-2}/t_{n-1})^2$")
ax.plot(ns, ru, "--", lw=2, label=r"rule  $(z^2+1)/[4(n-1)]$")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel("number of pairs  $n$")
ax.set_ylabel(r"break-even ICC  $\rho^*$")
ax.set_title(r"Pairing pays off only above $\rho^*$ (and $\rho^*\to 0$ slowly)")
ax.legend(frameon=False)
ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "threshold_vs_n.png"), dpi=150)
plt.close(fig)

# Fig 2: power & half-width curves with rho* marked
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
ax = axes[0]
ax.plot(rhos, pw_p, "-o", ms=3, label="paired")
ax.plot(rhos, pw_c, "-s", ms=3, label="complete rand.")
ax.axvline(rho_star_demo, color="k", ls=":", label=r"$\rho^*=%.3f$" % rho_star_demo)
ax.set_xlabel(r"ICC  $\rho$"); ax.set_ylabel("power")
ax.set_title(f"Power ($n={n_demo}$ pairs, $N={2*n_demo}$)")
ax.legend(frameon=False); ax.grid(alpha=0.25)
ax = axes[1]
ax.plot(rhos, hw_p, "-o", ms=3, label="paired")
ax.plot(rhos, hw_c, "-s", ms=3, label="complete rand.")
ax.axvline(rho_star_demo, color="k", ls=":", label=r"$\rho^*=%.3f$" % rho_star_demo)
ax.set_xlabel(r"ICC  $\rho$"); ax.set_ylabel("mean 95% CI half-width")
ax.set_title("Interval width")
ax.legend(frameon=False); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "power_width_crossing.png"), dpi=150)
plt.close(fig)

# Fig 3: rho*(b) non-monotonicity + Fisher contrast
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(bs_ok, rs_width, "-o", lw=2, label="width / power threshold")
ax.plot(bs_ok, rs_fisher, "-^", lw=2, label="Fisher / Cochran–Cox RE threshold")
ax.set_xlabel("block size  $b$  (finer blocking $\\leftarrow$)")
ax.set_ylabel(r"break-even ICC  $\rho^*$")
ax.set_title(r"Fixed $N=%d$: finer blocking demands higher ICC" % N_fixed)
ax.invert_xaxis()
ax.legend(frameon=False); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "threshold_vs_blocksize.png"), dpi=150)
plt.close(fig)

print("\nFigures written to", FIGDIR)
print("  threshold_vs_n.png, power_width_crossing.png, threshold_vs_blocksize.png")

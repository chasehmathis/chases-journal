"""
Averaging e-values cannot accumulate: the chi-square ceiling and a power reversal.

We work with the canonical Gaussian likelihood-ratio e-value. For a single study,
observe X and test H0: X ~ p = N(0,1) against the alternative q = N(delta, 1).
The likelihood-ratio e-variable is

    E = q(X)/p(X) = exp(delta * X - delta^2 / 2),

which is a valid e-value: E_{p}[E] = 1. Under the alternative (X ~ q) it has

    E_{q}[log E] = delta^2 / 2 = KL(q||p)      (per-study growth rate)
    E_{q}[E]     = exp(delta^2) = 1 + chi^2(q||p)   (alternative mean).

Given K i.i.d. studies E_1,...,E_K, the two canonical merged e-values are

    Product  P_K = prod_k E_k        (valid under independence)
    Average  A_K = (1/K) sum_k E_k   (valid under ARBITRARY dependence).

This script verifies:
  (1) the two constants (KL and 1+chi^2) for the Gaussian LR e-value;
  (2) the dichotomy: (1/K) log P_K -> KL, while A_K -> 1 + chi^2 (a finite ceiling);
  (3) validity of both merged e-values under the null (Type-I <= alpha);
  (4) the power reversal: the averaging test's power -> 0 for every level
      alpha < alpha* = 1/(1+chi^2), and -> 1 for alpha > alpha*, while the
      product test's power -> 1 at every fixed level.

All figures are written to figs/.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

SEED = 20260709
rng = np.random.default_rng(SEED)

FIGS = Path(__file__).resolve().parent / "figs"
FIGS.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Model constants
# ---------------------------------------------------------------------------
DELTA = 1.0                      # alternative shift; p=N(0,1), q=N(delta,1)
KL = DELTA**2 / 2.0              # per-study growth rate of the product
CHI2 = np.exp(DELTA**2) - 1.0    # chi-square divergence
CEIL = 1.0 + CHI2                # alternative mean of a single e-value = e^{delta^2}
ALPHA_STAR = 1.0 / CEIL          # phase-transition level = e^{-delta^2}

print("=" * 70)
print(f"Gaussian LR e-value:  p=N(0,1)  vs  q=N({DELTA},1)")
print(f"  KL(q||p)            = {KL:.6f}   (per-study product growth rate)")
print(f"  chi^2(q||p)         = {CHI2:.6f}")
print(f"  ceiling 1+chi^2     = {CEIL:.6f}   (= e^{{delta^2}} = {np.exp(DELTA**2):.6f})")
print(f"  alpha* = 1/(1+chi^2)= {ALPHA_STAR:.6f}")
print("=" * 70)


def draw_logE(n_rows, n_cols, mean):
    """logE_i for X ~ N(mean, 1):  logE = delta*X - delta^2/2."""
    X = rng.normal(mean, 1.0, size=(n_rows, n_cols))
    return DELTA * X - DELTA**2 / 2.0


# ---------------------------------------------------------------------------
# (1) Verify the two constants empirically
# ---------------------------------------------------------------------------
Nbig = 4_000_000
logE = draw_logE(1, Nbig, mean=DELTA).ravel()   # under the alternative
E = np.exp(logE)
print("\n[1] Monte-Carlo check of the two constants (alternative, "
      f"N={Nbig:,}):")
print(f"    mean(log E)  = {logE.mean():.5f}   (theory KL      = {KL:.5f})")
print(f"    mean(E)      = {E.mean():.5f}   (theory 1+chi^2 = {CEIL:.5f})")

# ---------------------------------------------------------------------------
# (2) Dichotomy: a single long path of P_K and A_K
# ---------------------------------------------------------------------------
Kmax_path = 500
logE_path = draw_logE(1, Kmax_path, mean=DELTA).ravel()
Ep = np.exp(logE_path)
ks = np.arange(1, Kmax_path + 1)
logP = np.cumsum(logE_path)                  # log of product
A = np.cumsum(Ep) / ks                        # running average

fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
ax[0].plot(ks, logP, lw=1.1, color="#1f77b4", label=r"$\log P_K$ (product)")
ax[0].plot(ks, ks * KL, "--", color="#111111", lw=1.2,
           label=r"$K\cdot \mathrm{KL}=K\delta^2/2$")
ax[0].set_xlabel("number of studies $K$")
ax[0].set_ylabel("log-evidence")
ax[0].set_title("Product: log-evidence grows linearly")
ax[0].legend(frameon=False, fontsize=9)

ax[1].plot(ks, A, lw=1.1, color="#d62728", label=r"$A_K$ (average)")
ax[1].axhline(CEIL, ls="--", color="#111111", lw=1.2,
              label=r"$1+\chi^2=e^{\delta^2}$")
ax[1].set_xlabel("number of studies $K$")
ax[1].set_ylabel("e-value $A_K$")
ax[1].set_title("Average: saturates at the $1+\\chi^2$ ceiling")
ax[1].legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(FIGS / "dichotomy.png", dpi=150)
plt.close(fig)
print("\n[2] Single-path dichotomy written to figs/dichotomy.png")
print(f"    log P_K at K=500 = {logP[-1]:.2f}   (theory ~ K*KL = {Kmax_path*KL:.2f})")
print(f"    A_K   at K=500   = {A[-1]:.4f}   (ceiling = {CEIL:.4f})")

# ---------------------------------------------------------------------------
# (3)+(4) Power / Type-I as a function of K, by Monte Carlo over reps
# ---------------------------------------------------------------------------
REPS = 20000
Kmax = 200
Kgrid = np.array([1, 2, 3, 5, 8, 12, 20, 30, 50, 75, 100, 150, 200])
alphas = [0.5, 0.4, 0.3, 0.1, 0.05]   # straddle alpha* = 1/e ~ 0.3679


def power_curves(mean):
    """For X~N(mean,1): return dict alpha -> (avg_power[K], prod_power[K]) over Kgrid.
    Power = P(merged e-value >= 1/alpha)."""
    logEm = draw_logE(REPS, Kmax, mean=mean)     # reps x Kmax
    Em = np.exp(logEm)
    cum_logP = np.cumsum(logEm, axis=1)          # log product, reps x Kmax
    cum_A = np.cumsum(Em, axis=1) / np.arange(1, Kmax + 1)  # average, reps x Kmax
    out = {}
    for a in alphas:
        thr = np.log(1.0 / a)
        prodP = (cum_logP >= thr).mean(axis=0)   # per-K power of product test
        avgP = (cum_A >= 1.0 / a).mean(axis=0)    # per-K power of average test
        out[a] = (avgP[Kgrid - 1], prodP[Kgrid - 1])
    return out


alt = power_curves(mean=DELTA)   # power (alternative true)
nul = power_curves(mean=0.0)     # Type-I (null true)

print(f"\n[3] Type-I error under the null (X~N(0,1)), REPS={REPS:,}:")
print("     alpha   avg P(A_K>=1/a) [max over K]   prod P(P_K>=1/a) [max over K]")
for a in alphas:
    avgP, prodP = nul[a]
    print(f"     {a:<6}  {avgP.max():.4f}"
          f"                        {prodP.max():.4f}   (<= {a})")

print(f"\n[4] Power under the alternative (X~N({DELTA},1)); alpha* = {ALPHA_STAR:.4f}")
print("     Averaging test power P(A_K >= 1/alpha):")
print("        K:   " + "  ".join(f"{k:>5d}" for k in Kgrid))
for a in alphas:
    tag = "below a*" if a < ALPHA_STAR else "above a*"
    avgP, _ = alt[a]
    print(f"     a={a:<4} ({tag}): " + "  ".join(f"{v:5.3f}" for v in avgP))
print("     Product test power P(P_K >= 1/alpha):")
for a in alphas:
    _, prodP = alt[a]
    print(f"     a={a:<4}         : " + "  ".join(f"{v:5.3f}" for v in prodP))

# ---- Figure: power vs K, averaging vs product -----------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9), sharey=True)
cmap = plt.cm.viridis(np.linspace(0.1, 0.85, len(alphas)))
for a, c in zip(alphas, cmap):
    avgP, prodP = alt[a]
    style = "-o" if a < ALPHA_STAR else "-s"
    ax[0].plot(Kgrid, avgP, style, color=c, ms=4, lw=1.2,
               label=rf"$\alpha={a}$" + (" $<\\alpha^*$" if a < ALPHA_STAR else " $>\\alpha^*$"))
    ax[1].plot(Kgrid, prodP, "-o", color=c, ms=4, lw=1.2, label=rf"$\alpha={a}$")
ax[0].set_title(r"Averaging test: power reverses at $\alpha^*=1/(1+\chi^2)$")
ax[0].set_xlabel("number of studies $K$")
ax[0].set_ylabel(r"power $=\Pr(\text{reject})$")
ax[0].set_xscale("log")
ax[0].legend(frameon=False, fontsize=8, loc="center right")
ax[1].set_title("Product test: power $\\to 1$ at every level")
ax[1].set_xlabel("number of studies $K$")
ax[1].set_xscale("log")
ax[1].legend(frameon=False, fontsize=8, loc="lower right")
for a_ in ax:
    a_.set_ylim(-0.03, 1.03)
fig.tight_layout()
fig.savefig(FIGS / "power_reversal.png", dpi=150)
plt.close(fig)
print("\n    Power figure written to figs/power_reversal.png")

# ---------------------------------------------------------------------------
# (5) The ceiling as a function of delta: A_infty = e^{delta^2}
# ---------------------------------------------------------------------------
deltas = np.linspace(0.25, 1.75, 13)
Kbig = 4000
emp_ceiling = []
for d in deltas:
    X = rng.normal(d, 1.0, size=(400, Kbig))
    Ed = np.exp(d * X - d**2 / 2.0)
    emp_ceiling.append((np.cumsum(Ed, axis=1)[:, -1] / Kbig).mean())
emp_ceiling = np.array(emp_ceiling)
theory_ceiling = np.exp(deltas**2)

fig, ax = plt.subplots(figsize=(5.4, 3.9))
ax.plot(deltas, theory_ceiling, "--", color="#111111", lw=1.3,
        label=r"theory $1+\chi^2=e^{\delta^2}$")
ax.plot(deltas, emp_ceiling, "o", color="#d62728", ms=5,
        label=rf"$A_{{{Kbig}}}$ (Monte Carlo)")
ax.set_yscale("log")
ax.set_xlabel(r"alternative shift $\delta$")
ax.set_ylabel("averaging ceiling")
ax.set_title(r"The average saturates at $e^{\delta^2}$, no matter how large $K$")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(FIGS / "ceiling_vs_delta.png", dpi=150)
plt.close(fig)
print("[5] Ceiling-vs-delta written to figs/ceiling_vs_delta.png")
print("    delta :", "  ".join(f"{d:.2f}" for d in deltas))
print("    A_big :", "  ".join(f"{v:6.2f}" for v in emp_ceiling))
print("    e^d^2 :", "  ".join(f"{v:6.2f}" for v in theory_ceiling))
print("\nDone.")

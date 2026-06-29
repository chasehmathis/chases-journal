"""
Effective dimension of equicorrelated multiple testing.

Companion code to the note. Two computations:

  (1) EXACT FWER of one-sided Bonferroni for m equicorrelated standard-normal
      test statistics under the global null, via a 1-D Gaussian quadrature
      (conditioning on the common factor Z). No Monte Carlo error.
  (2) A Monte Carlo check of that quadrature at moderate m, with error bars.

We then verify, against these exact numbers:

  Result 1 (overconservativeness).  At the nominal Bonferroni threshold
  t_m = Phibar^{-1}(alpha/m),
        FWER(m, rho)  =  m^{-beta(rho) + o(1)},     beta(rho) = (1 - sqrt(1-rho))^2 / rho.
  We estimate the local log-log slope of FWER(m) and compare to -beta(rho).

  Result 2 (effective multiplicity).  The number m_eff = alpha / Phibar(t*),
  where t* solves FWER(m, rho; t*) = alpha, satisfies
        log m_eff = (1 - rho) log m + O(sqrt(log m)),
  i.e. the *leading* exponent is 1 - rho (convergence is slow, O(1/sqrt(log m))).

Figures written to figs/.  Seed fixed for the Monte Carlo.
"""

import numpy as np
from scipy.stats import norm
from scipy.integrate import quad
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG = np.random.default_rng(20260629)
ALPHA = 0.05


# ----------------------------------------------------------------------------
# Exact equicorrelated FWER by 1-D quadrature.
#   X_i = sqrt(rho) Z + sqrt(1-rho) eps_i,  Z, eps_i iid N(0,1).
#   Given Z=z, the X_i are iid with P(X_i <= t | z) = Phi((t - sqrt(rho) z)/sqrt(1-rho)).
#   FWER = P(max X_i > t) = 1 - E_Z[ Phi((t - sqrt(rho)Z)/sqrt(1-rho))^m ].
# ----------------------------------------------------------------------------
def fwer_exact(m, rho, t):
    if rho <= 0.0:
        return 1.0 - norm.cdf(t) ** m
    a, b = np.sqrt(rho), np.sqrt(1.0 - rho)

    def integrand(z):
        u = (t - a * z) / b
        return norm.pdf(z) * np.exp(m * norm.logcdf(u))

    val, _ = quad(integrand, -12.0, 16.0, limit=800)
    return 1.0 - val


def bonferroni_threshold(m, alpha=ALPHA):
    return norm.isf(alpha / m)


def beta_exponent(rho):
    """Predicted FWER-decay exponent."""
    return (1.0 - np.sqrt(1.0 - rho)) ** 2 / rho


# ----------------------------------------------------------------------------
# Monte Carlo FWER (validation of the quadrature).
# ----------------------------------------------------------------------------
def fwer_mc(m, rho, t, reps, max_cells=10_000_000):
    a, b = np.sqrt(rho), np.sqrt(1.0 - rho)
    batch = max(1, max_cells // m)  # cap memory: batch * m cells per draw
    hits = 0
    done = 0
    while done < reps:
        n = min(batch, reps - done)
        Z = RNG.standard_normal((n, 1))
        E = RNG.standard_normal((n, m))
        X = a * Z + b * E
        hits += int(np.sum(X.max(axis=1) > t))
        done += n
    p = hits / reps
    se = np.sqrt(p * (1 - p) / reps)
    return p, se


# ============================================================================
# Result 1: FWER decay exponent.
# ============================================================================
print("=" * 74)
print("Result 1: nominal Bonferroni FWER decays as m^{-beta(rho)}")
print("=" * 74)
rhos1 = [0.1, 0.3, 0.5, 0.8]
ms = np.array([10 ** k for k in range(2, 9)])  # 1e2 ... 1e8

fwer_curves = {}
for rho in rhos1:
    fw = np.array([fwer_exact(int(m), rho, bonferroni_threshold(int(m))) for m in ms])
    fwer_curves[rho] = fw
    local = (np.log(fw[-1]) - np.log(fw[-2])) / (np.log(ms[-1]) - np.log(ms[-2]))
    overall = np.polyfit(np.log(ms), np.log(fw), 1)[0]
    print(
        f"  rho={rho:>4}:  beta(rho)={beta_exponent(rho):.4f}   "
        f"-(local slope 1e7->1e8)={-local:.4f}   "
        f"FWER@1e6={fw[4]:.3e}   FWER@1e8={fw[-1]:.3e}"
    )

# Monte Carlo validation at moderate m.
print("\nMonte Carlo validation of the quadrature (FWER, global null):")
mc_rows = []
for (m, rho, reps) in [(1000, 0.3, 400000), (10000, 0.5, 400000), (50000, 0.8, 200000)]:
    t = bonferroni_threshold(m)
    exact = fwer_exact(m, rho, t)
    p, se = fwer_mc(m, rho, t, reps=reps)
    mc_rows.append((m, rho, exact, p, se))
    print(
        f"  m={m:>6}, rho={rho}:  quadrature={exact:.4e}   "
        f"MC={p:.4e} +/- {1.96*se:.1e} (95% CI)"
    )

# ============================================================================
# Result 2: effective multiplicity to recover FWER = alpha.
# ============================================================================
print("\n" + "=" * 74)
print("Result 2: effective number of tests m_eff ~ m^{1-rho} (slow convergence)")
print("=" * 74)
rhos2 = [0.2, 0.5, 0.8]
ms2 = np.array([10 ** k for k in range(3, 10)])  # 1e3 ... 1e9
meff_curves = {}
for rho in rhos2:
    meff = []
    for m in ms2:
        m = int(m)
        f = lambda t: fwer_exact(m, rho, t) - ALPHA
        tstar = brentq(f, 0.5, bonferroni_threshold(m) + 3.0, xtol=1e-7)
        meff.append(ALPHA / norm.sf(tstar))
    meff = np.array(meff)
    meff_curves[rho] = meff
    local = (np.log(meff[-1]) - np.log(meff[-2])) / (np.log(ms2[-1]) - np.log(ms2[-2]))
    print(
        f"  rho={rho}:  predicted leading exponent 1-rho={1-rho:.3f}   "
        f"local slope(1e8->1e9)={local:.4f}   m_eff@1e6={meff[3]:.3e}  (of m=1e6)"
    )

# ============================================================================
# Figures
# ============================================================================
# Figure 1: FWER decay, log-log, exact curves + theory slopes + MC points.
fig, ax = plt.subplots(figsize=(7.2, 5.0))
colors = plt.cm.viridis(np.linspace(0.05, 0.8, len(rhos1)))
for rho, c in zip(rhos1, colors):
    fw = fwer_curves[rho]
    ax.loglog(ms, fw, "o-", color=c, lw=1.8, ms=4, label=f"$\\rho={rho}$ (exact)")
    # theory slope line anchored at m=1e4
    anchor_i = 2
    ref = fw[anchor_i] * (ms / ms[anchor_i]) ** (-beta_exponent(rho))
    ax.loglog(ms, ref, "--", color=c, lw=1.0, alpha=0.7)
for (m, rho, exact, p, se) in mc_rows:
    ax.errorbar(m, p, yerr=1.96 * se, fmt="s", color="crimson", ms=6,
                capsize=3, zorder=5)
ax.axhline(ALPHA, color="grey", ls=":", lw=1.0)
ax.text(ms[0], ALPHA * 1.15, r"nominal $\alpha=0.05$", color="grey", fontsize=9)
ax.errorbar([], [], yerr=[], fmt="s", color="crimson", capsize=3,
            label="Monte Carlo (95% CI)")
ax.set_xlabel("number of tests $m$")
ax.set_ylabel(r"FWER at nominal Bonferroni threshold $\alpha/m$")
ax.set_title("Nominal Bonferroni is polynomially overconservative under equicorrelation\n"
             r"dashed: predicted slope $-\beta(\rho)=-(1-\sqrt{1-\rho})^2/\rho$")
ax.legend(fontsize=8.5, loc="lower left")
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig("figs/fwer_decay.png", dpi=140)
plt.close(fig)

# Figure 2: two panels — (a) exponent vs rho; (b) effective multiplicity.
fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.7))

# (a) beta(rho) curve and the effective-multiplicity leading exponent 1-rho.
rr = np.linspace(0.001, 0.999, 400)
axL.plot(rr, beta_exponent(rr), "b-", lw=2.2,
         label=r"$\beta(\rho)=\dfrac{(1-\sqrt{1-\rho})^2}{\rho}$  (FWER decay)")
axL.plot(rr, 1 - rr, "g-", lw=2.2,
         label=r"$1-\rho$  (effective multiplicity exponent)")
# overlay empirically estimated FWER-decay exponents (local slopes) at finite m
emp_rhos = [0.1, 0.3, 0.5, 0.8]
emp_beta = []
for rho in emp_rhos:
    fw = fwer_curves[rho]
    s = -(np.log(fw[-1]) - np.log(fw[-2])) / (np.log(ms[-1]) - np.log(ms[-2]))
    emp_beta.append(s)
axL.plot(emp_rhos, emp_beta, "bo", ms=7, mfc="white", mec="b",
         label=r"empirical $\beta$ (slope $10^7\!\to\!10^8$)")
axL.set_xlabel(r"equicorrelation $\rho$")
axL.set_ylabel("exponent")
axL.set_title("Two sharp exponents of the collapse")
axL.legend(fontsize=9, loc="upper center")
axL.grid(True, alpha=0.25)
axL.set_xlim(0, 1)
axL.set_ylim(0, 1)

# (b) effective multiplicity vs m with m^{1-rho} reference.
colors2 = plt.cm.plasma(np.linspace(0.1, 0.75, len(rhos2)))
for rho, c in zip(rhos2, colors2):
    me = meff_curves[rho]
    axR.loglog(ms2, me, "o-", color=c, lw=1.8, ms=4, label=f"$\\rho={rho}$: $m_{{eff}}$")
    ref = me[2] * (ms2 / ms2[2]) ** (1 - rho)
    axR.loglog(ms2, ref, "--", color=c, lw=1.0, alpha=0.7)
axR.loglog(ms2, ms2, "k:", lw=1.0, label=r"$m_{eff}=m$ (independence)")
axR.set_xlabel("number of tests $m$")
axR.set_ylabel(r"effective multiplicity $m_{eff}$ (for FWER $=\alpha$)")
axR.set_title(r"$m_{eff}\sim m^{1-\rho}$ (dashed); grows with $m$, not constant")
axR.legend(fontsize=8.5, loc="upper left")
axR.grid(True, which="both", alpha=0.25)

fig.tight_layout()
fig.savefig("figs/exponents.png", dpi=140)
plt.close(fig)

print("\nFigures written: figs/fwer_decay.png, figs/exponents.png")

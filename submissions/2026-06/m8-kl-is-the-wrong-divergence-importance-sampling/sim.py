"""
Simulations for "KL is the wrong divergence for importance sampling:
a sharp factor-of-two for Gaussians."

Setup. Proposal q = N(0,1). Target p = N(0, s2) (s2 = sigma^2), scale family.
Self-normalized importance sampling (SNIS) estimates a target expectation
E_p[f] using x_i ~ q and weights w_i = p(x_i)/q(x_i):

        hat = sum_i w_i f(x_i) / sum_i w_i.

We estimate the second moment, f(x) = x^2, whose true value under p is s2.

Theory we check:
  * Closed form  1 + chi^2(p||q) = E_q[w^2] = 1 / (sigma * sqrt(2 - sigma^2)),
    finite iff sigma^2 < 2; the importance weights have INFINITE variance for
    sigma^2 >= 2.
  * KL(p||q) = (sigma^2 - 1 - ln sigma^2)/2 stays small (~0.153 at sigma^2=2),
    so a KL-based sample-size rule (n ~ exp(KL)) is blind to the blow-up.
  * Consequences: for sigma^2 < 2 the SNIS error decays like 1/sqrt(n);
    for sigma^2 >= 2 it does not, and the Kish effective sample size (ESS)
    gives false confidence.

All figures saved to figs/. Seed fixed for reproducibility.
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG = np.random.default_rng(20260625)
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)


# ----------------------------------------------------------------------
# Closed-form quantities (scale family: q=N(0,1), p=N(0,s2))
# ----------------------------------------------------------------------
def one_plus_chi2_scale(s2):
    """E_q[w^2] = 1/(sigma*sqrt(2-s2)) for s2<2, else +inf."""
    s2 = np.asarray(s2, dtype=float)
    out = np.full_like(s2, np.inf)
    mask = s2 < 2.0
    sig = np.sqrt(s2[mask])
    out[mask] = 1.0 / (sig * np.sqrt(2.0 - s2[mask]))
    return out


def kl_scale(s2):
    """KL(N(0,s2) || N(0,1)) = (s2 - 1 - ln s2)/2."""
    s2 = np.asarray(s2, dtype=float)
    return 0.5 * (s2 - 1.0 - np.log(s2))


def log_weights(x, s2):
    """log p(x)/q(x) for p=N(0,s2), q=N(0,1)."""
    return -0.5 * np.log(s2) - 0.5 * x * x * (1.0 / s2 - 1.0)


# ----------------------------------------------------------------------
# (0) Numerically validate the closed form for 1 + chi^2
# ----------------------------------------------------------------------
def validate_closed_form():
    print("=== Validation of 1 + chi^2 = E_q[w^2] = 1/(sigma*sqrt(2-s2)) ===")
    N = 20_000_000
    x = RNG.standard_normal(N)
    for s2 in [0.5, 1.0, 1.5, 1.8]:
        lw = log_weights(x, s2)
        mc = np.mean(np.exp(2.0 * lw))  # E_q[w^2]
        cf = float(one_plus_chi2_scale(np.array([s2]))[0])
        print(f"  s2={s2:>4}:  MC E_q[w^2]={mc:10.4f}   closed form={cf:10.4f}")
    print()


# ----------------------------------------------------------------------
# (1) Divergence figure: chi^2 blows up at s2=2, KL stays tiny
# ----------------------------------------------------------------------
def fig_divergence():
    s2 = np.linspace(0.2, 3.0, 600)
    chi2 = one_plus_chi2_scale(s2) - 1.0
    kl = kl_scale(s2)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    ax[0].plot(s2, kl, color="C0", lw=2)
    ax[0].axvline(2.0, color="k", ls="--", lw=1)
    ax[0].set_title("KL$(p\\,\\|\\,q)$  (bounded, smooth through $\\sigma^2=2$)")
    ax[0].set_xlabel("$\\sigma^2$ (target variance)")
    ax[0].set_ylabel("KL divergence (nats)")
    ax[0].annotate("KL $\\approx 0.153$\nat $\\sigma^2=2$", xy=(2.0, kl_scale(np.array([2.0]))[0]),
                   xytext=(2.1, 0.6), arrowprops=dict(arrowstyle="->"))
    ax[0].set_ylim(0, 1.2)

    chi2_plot = np.where(s2 < 2.0, chi2, np.nan)
    ax[1].plot(s2, chi2_plot, color="C3", lw=2)
    ax[1].axvline(2.0, color="k", ls="--", lw=1)
    ax[1].set_title("$\\chi^2(p\\,\\|\\,q)=E_q[w^2]-1$  (diverges at $\\sigma^2=2$)")
    ax[1].set_xlabel("$\\sigma^2$ (target variance)")
    ax[1].set_ylabel("$\\chi^2$ divergence")
    ax[1].set_ylim(0, 8)
    ax[1].annotate("$\\to\\infty$", xy=(1.95, 7), fontsize=13, color="C3")

    fig.suptitle("Same proposal, same target family: KL is blind to the variance blow-up")
    fig.tight_layout()
    p = os.path.join(FIGS, "divergence.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"wrote {p}")


# ----------------------------------------------------------------------
# (2) SNIS estimate of E_p[X^2]=s2: error vs n, finite vs infinite variance
# ----------------------------------------------------------------------
def snis_estimate(s2, n, reps):
    """Return array (reps,) of SNIS estimates of E_p[X^2] with sample size n."""
    est = np.empty(reps)
    for r in range(reps):
        x = RNG.standard_normal(n)
        lw = log_weights(x, s2)
        lw -= lw.max()  # stabilize
        w = np.exp(lw)
        est[r] = np.sum(w * x * x) / np.sum(w)
    return est


def fig_rmse_vs_n():
    ns = np.array([100, 300, 1000, 3000, 10000, 30000, 100000])
    reps = 400
    scenarios = {1.5: "C0", 2.5: "C3", 3.0: "C1"}
    print("=== SNIS RMSE vs n (target E_p[X^2]=s2) ===")
    fig, ax = plt.subplots(figsize=(7, 5))
    for s2, color in scenarios.items():
        rmse = []
        for n in ns:
            est = snis_estimate(s2, int(n), reps)
            rmse.append(np.sqrt(np.mean((est - s2) ** 2)))
        rmse = np.array(rmse)
        finite = "finite var" if s2 < 2 else "INFINITE var"
        ax.loglog(ns, rmse, "o-", color=color,
                  label=f"$\\sigma^2={s2}$ ({finite})")
        print(f"  s2={s2}: RMSE @ n=100 -> {rmse[0]:.3f},  @ n=100000 -> {rmse[-1]:.3f}"
              f"   ratio={rmse[0]/rmse[-1]:.1f} (1/sqrt(1000)=31.6 ideal)")
    # reference 1/sqrt(n) slope
    ref = 1.5 * (ns / ns[0]) ** -0.5 * 0.5
    ax.loglog(ns, ref, "k--", lw=1, label="$1/\\sqrt{n}$ reference")
    ax.set_xlabel("sample size $n$")
    ax.set_ylabel("RMSE of SNIS estimate of $E_p[X^2]$")
    ax.set_title("Below $\\sigma^2=2$: $1/\\sqrt{n}$. At/above: error refuses to shrink.")
    ax.legend()
    fig.tight_layout()
    p = os.path.join(FIGS, "rmse_vs_n.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"wrote {p}\n")


# ----------------------------------------------------------------------
# (3) Running estimate: a single long run, finite vs infinite variance
# ----------------------------------------------------------------------
def fig_running_estimate():
    n = 200_000
    fig, ax = plt.subplots(figsize=(7, 5))
    grid = np.unique(np.logspace(1, np.log10(n), 300).astype(int))
    for s2, color in [(1.5, "C0"), (2.5, "C3"), (3.0, "C1")]:
        x = RNG.standard_normal(n)
        lw = log_weights(x, s2)
        lw -= lw.max()
        w = np.exp(lw)
        num = np.cumsum(w * x * x)
        den = np.cumsum(w)
        run = num[grid - 1] / den[grid - 1]
        finite = "finite var" if s2 < 2 else "infinite var"
        ax.semilogx(grid, run, color=color, lw=1.4,
                    label=f"$\\sigma^2={s2}$ ({finite}), target$={s2}$")
        ax.axhline(s2, color=color, ls=":", lw=1)
    ax.set_xlabel("sample size $n$ (running SNIS)")
    ax.set_ylabel("running estimate of $E_p[X^2]$")
    ax.set_title("One long run: smooth convergence vs. jump-driven wandering")
    ax.legend()
    fig.tight_layout()
    p = os.path.join(FIGS, "running_estimate.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"wrote {p}")


# ----------------------------------------------------------------------
# (4) Kish ESS gives false confidence in the infinite-variance regime
# ----------------------------------------------------------------------
def fig_ess_distribution():
    n = 5000
    reps = 3000
    print("=== Kish ESS distribution (n=5000), relative ESS = ESS/n ===")
    fig, ax = plt.subplots(figsize=(7, 5))
    for s2, color in [(1.5, "C0"), (3.0, "C3")]:
        ress = np.empty(reps)
        for r in range(reps):
            x = RNG.standard_normal(n)
            lw = log_weights(x, s2)
            lw -= lw.max()
            w = np.exp(lw)
            ress[r] = (w.sum() ** 2) / (n * np.sum(w * w))
        finite = "finite var" if s2 < 2 else "infinite var (true rel-ESS $\\to 0$)"
        ax.hist(ress, bins=60, range=(0, 1), alpha=0.6, color=color,
                density=True, label=f"$\\sigma^2={s2}$ ({finite})")
        print(f"  s2={s2}: mean rel-ESS={ress.mean():.3f}, "
              f"median={np.median(ress):.3f}, "
              f"P(rel-ESS>0.3)={np.mean(ress>0.3):.3f}")
    ax.set_xlabel("Kish relative ESS  $=(\\sum w)^2 / (n\\sum w^2)$")
    ax.set_ylabel("density over independent runs")
    ax.set_title("The ESS diagnostic looks healthy even when variance is infinite")
    ax.legend()
    fig.tight_layout()
    p = os.path.join(FIGS, "ess_distribution.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print(f"wrote {p}\n")


# ----------------------------------------------------------------------
# (5) Location vs scale dichotomy: report the key numbers
# ----------------------------------------------------------------------
def report_dichotomy():
    print("=== Location vs scale: KL and log(1+chi^2) ===")
    print("  SCALE family q=N(0,1), p=N(0,s2):")
    for s2 in [1.5, 1.9, 1.99]:
        kl = float(kl_scale(np.array([s2]))[0])
        lc = np.log(float(one_plus_chi2_scale(np.array([s2]))[0]))
        print(f"    s2={s2:>5}: KL={kl:.4f},  log(1+chi^2)={lc:.4f}")
    print("    s2>=2   : KL finite (<=0.153 region near 2), log(1+chi^2)=+inf")
    print("  LOCATION family q=N(0,1), p=N(mu,1):  KL=mu^2/2, log(1+chi^2)=mu^2 = 2*KL")
    for mu in [1.0, 2.0, 3.0]:
        print(f"    mu={mu}: KL={mu**2/2:.3f},  log(1+chi^2)={mu**2:.3f}")
    print()


if __name__ == "__main__":
    validate_closed_form()
    report_dichotomy()
    fig_divergence()
    fig_rmse_vs_n()
    fig_running_estimate()
    fig_ess_distribution()
    print("done.")

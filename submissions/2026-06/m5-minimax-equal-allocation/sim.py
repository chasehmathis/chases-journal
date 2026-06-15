"""
Simulations for "Fifty-fifty is not naive: equal allocation is the minimax
design for treatment-effect estimation under unknown variances."

Two arms, outcome variances (sigma0^2, sigma1^2). Allocate fractions (1-p, p)
of n units to control/treatment. The difference-in-means estimator of the ATE
has variance

    V(p) = sigma0^2/(n(1-p)) + sigma1^2/(n p).

Write t = sigma1/sigma0 (the variance-RATIO, in sd units) and normalize
sigma0 = 1. Then (dropping the 1/n) Vbar(p;t) = 1/(1-p) + t^2/p, the Neyman
optimum is p*(t) = t/(1+t) with Vbar(p*) = (1+t)^2, and the relative
efficiency of design p against the oracle is

    rho(p,t) = Vbar(p*;t) / Vbar(p;t) = (1+t)^2 / ( 1/(1-p) + t^2/p )  in (0,1].

This script:
  (A) plots the worst-case (over t) relative efficiency g(p)=inf_t rho(p,t)
      and confirms the maximin design is p=1/2 with guarantee exactly 1/2;
  (B) Monte-Carlo checks the factor-of-two and the analytic V(p);
  (C) studies a two-stage adaptive-Neyman design: does estimating the
      variances recover the oracle, and when does adaptivity backfire?
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260615)
FIGS = __import__("os").path.join(__import__("os").path.dirname(__file__), "figs")
__import__("os").makedirs(FIGS, exist_ok=True)


# ----------------------------------------------------------------------
# Analytic helpers (sigma0 = 1, sigma1 = t)
# ----------------------------------------------------------------------
def Vbar(p, t):
    return 1.0 / (1.0 - p) + t**2 / p

def rho(p, t):
    return (1.0 + t)**2 / Vbar(p, t)

def pstar(t):
    return t / (1.0 + t)


# ----------------------------------------------------------------------
# (A) Worst-case relative efficiency g(p) = inf_t rho(p,t)
# ----------------------------------------------------------------------
def worst_case_efficiency():
    p_grid = np.linspace(0.02, 0.98, 193)
    t_grid = np.concatenate([np.logspace(-4, 0, 400), np.logspace(0, 4, 400)])
    g = np.array([rho(p, t_grid).min() for p in p_grid])

    p_best = p_grid[np.argmax(g)]
    print("=== (A) Minimax (maximin) design ===")
    print(f"  argmax_p inf_t rho(p,t)  = {p_best:.4f}  (theory: 0.5)")
    print(f"  guaranteed efficiency    = {g.max():.4f}  (theory: 0.5)")
    # compare to closed form min(p,1-p)
    closed = np.minimum(p_grid, 1 - p_grid)
    print(f"  max |g(p) - min(p,1-p)|  = {np.abs(g - closed).max():.4f}")

    # worst-case efficiency of a few named designs
    for p in (0.5, 0.6, 0.7, 0.8):
        gp = rho(p, t_grid).min()
        print(f"  inf_t rho(p={p:.2f}) = {gp:.4f}   (min(p,1-p)={min(p,1-p):.3f})")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.0))
    ax1.plot(p_grid, g, lw=2.2, color="C0", label=r"$g(p)=\inf_t\,\rho(p,t)$")
    ax1.plot(p_grid, closed, "--", color="C3", lw=1.4,
             label=r"$\min(p,1-p)$")
    ax1.axvline(0.5, color="0.6", lw=1, ls=":")
    ax1.scatter([0.5], [0.5], color="C0", zorder=5)
    ax1.annotate("minimax design\n(p=1/2, guarantee 1/2)", (0.5, 0.5),
                 xytext=(0.13, 0.30),
                 arrowprops=dict(arrowstyle="->", color="0.4"), fontsize=9)
    ax1.set_xlabel("treatment fraction $p$")
    ax1.set_ylabel("worst-case relative efficiency")
    ax1.set_ylim(0, 1.02)
    ax1.legend(loc="upper right", fontsize=9)
    ax1.set_title("(a) Equal allocation is maximin")

    # rho(p,t) vs t for several designs
    tt = np.logspace(-1.3, 1.3, 300)
    for p, c in zip((0.5, 0.7, 0.9), ("C0", "C1", "C2")):
        ax2.plot(tt, rho(p, tt), color=c, lw=2, label=f"$p={p}$")
    ax2.axhline(0.5, color="0.6", lw=1, ls=":")
    ax2.set_xscale("log")
    ax2.set_xlabel(r"variance ratio $t=\sigma_1/\sigma_0$")
    ax2.set_ylabel(r"relative efficiency $\rho(p,t)$")
    ax2.set_ylim(0, 1.02)
    ax2.legend(loc="lower center", fontsize=9)
    ax2.set_title("(b) Efficiency profiles")
    fig.tight_layout()
    fig.savefig(f"{FIGS}/efficiency.png", dpi=140)
    plt.close(fig)
    return p_best, g.max()


# ----------------------------------------------------------------------
# (B) Monte-Carlo check of V(p) and the factor of two
# ----------------------------------------------------------------------
def mc_check_factor_two(n=2000, R=40000):
    print("\n=== (B) Monte-Carlo check of V(p) (n=%d, R=%d) ===" % (n, R))
    s0, s1 = 1.0, 4.0          # extreme ratio t=4
    t = s1 / s0
    for label, p in [("equal p=1/2", 0.5), ("Neyman p*", pstar(t))]:
        n1 = int(round(p * n)); n0 = n - n1
        est = np.empty(R)
        for r in range(R):
            y0 = rng.normal(0.0, s0, n0)
            y1 = rng.normal(0.0, s1, n1)
            est[r] = y1.mean() - y0.mean()
        emp = est.var()
        ana = s0**2 / n0 + s1**2 / n1
        print(f"  {label:12s} n1={n1:4d}: empVar={emp:.5f}  analytic={ana:.5f}")
    R_ratio = 2 * (s0**2 + s1**2) / (s0 + s1)**2
    print(f"  V(1/2)/V(p*) analytic ratio = {R_ratio:.4f}  (sup is 2)")


# ----------------------------------------------------------------------
# (C) Two-stage adaptive Neyman design
# ----------------------------------------------------------------------
def two_stage_once(n, s0, s1, f):
    """One run; returns the difference-in-means ATE estimate (true ATE = 0)."""
    n_p = max(4, int(round(f * n)))
    n_p += n_p % 2                      # even, split equally
    h = n_p // 2
    y0p = rng.normal(0.0, s0, h)
    y1p = rng.normal(0.0, s1, h)
    sd0 = max(y0p.std(ddof=1), 1e-6)
    sd1 = max(y1p.std(ddof=1), 1e-6)
    phat = sd1 / (sd0 + sd1)
    n_r = n - n_p
    n1r = int(round(phat * n_r)); n0r = n_r - n1r
    y0r = rng.normal(0.0, s0, n0r)
    y1r = rng.normal(0.0, s1, n1r)
    y0 = np.concatenate([y0p, y0r]); y1 = np.concatenate([y1p, y1r])
    return y1.mean() - y0.mean()


def fixed_once(n, s0, s1, p):
    n1 = int(round(p * n)); n0 = n - n1
    y0 = rng.normal(0.0, s0, n0); y1 = rng.normal(0.0, s1, n1)
    return y1.mean() - y0.mean()


def adaptive_study():
    print("\n=== (C) Two-stage adaptive Neyman ===")
    R = 20000
    f = 0.30                                   # pilot fraction
    ns = [40, 80, 160, 320, 640, 1280]
    results = {}
    for t in (1.0, 3.0, 9.0):
        s0, s1 = 1.0, t
        oracle = (s0 + s1)**2                   # n * Vbar(p*)
        rows = []
        for n in ns:
            v_eq = np.var([fixed_once(n, s0, s1, 0.5) for _ in range(R)]) * n
            v_or = np.var([fixed_once(n, s0, s1, pstar(t)) for _ in range(R)]) * n
            v_ad = np.var([two_stage_once(n, s0, s1, f) for _ in range(R)]) * n
            # relative efficiency = oracle variance / achieved variance
            rows.append((n, oracle / v_eq, oracle / v_or, oracle / v_ad))
            print(f"  t={t:>4}  n={n:>5}  effEqual={oracle/v_eq:.3f} "
                  f" effOracleMC={oracle/v_or:.3f}  effAdaptive={oracle/v_ad:.3f}")
        results[t] = np.array(rows)
    print(f"  (analytic equal-allocation efficiency: t=1 -> "
          f"{rho(0.5,1.0):.3f}, t=3 -> {rho(0.5,3.0):.3f}, "
          f"t=9 -> {rho(0.5,9.0):.3f})")

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0), sharey=True)
    for ax, t in zip(axes, (1.0, 3.0, 9.0)):
        rr = results[t]
        ax.axhline(1.0, color="0.7", lw=1, ls="--")
        ax.axhline(rho(0.5, t), color="C0", lw=1, ls=":")
        ax.plot(rr[:, 0], rr[:, 1], "o-", color="C0", label="equal (1/2)")
        ax.plot(rr[:, 0], rr[:, 3], "s-", color="C1", label="two-stage adaptive")
        ax.plot(rr[:, 0], rr[:, 2], "^--", color="0.5", label="oracle Neyman (MC)")
        ax.set_xscale("log")
        ax.set_xlabel("total sample size $n$")
        ax.set_title(f"$t=\\sigma_1/\\sigma_0={t:.0f}$  ($p^*={pstar(t):.2f}$)")
        ax.legend(loc="lower right", fontsize=8.5)
    axes[0].set_ylabel("efficiency vs oracle Neyman")
    fig.suptitle(f"Two-stage adaptive Neyman (pilot fraction f={f}), "
                 f"R={R} reps", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{FIGS}/adaptive.png", dpi=140)
    plt.close(fig)
    return results


if __name__ == "__main__":
    worst_case_efficiency()
    mc_check_factor_two()
    adaptive_study()
    print("\nDone. Figures in figs/.")

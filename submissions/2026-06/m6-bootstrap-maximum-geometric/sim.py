"""
The bootstrap of the sample maximum: an explicit compound-geometric limit law.

Data: X_1,...,X_n iid Uniform[0, theta], theta = 1.  M_n = max X_i.
True pivot: n(theta - M_n)/theta -> Exp(1).

We verify three things and produce two figures:

  (A) RANK deficit.  Let M_n* = max of a nonparametric bootstrap resample, and
      D = (rank of M_n among the order stats) - (rank of M_n*).  We claim the
      bootstrap law of D converges to Geometric(1 - 1/e) on {0,1,2,...}:
          P(D = j) -> (1 - e^{-1}) e^{-j}.
      In particular P(D = 0) = P(M_n* = M_n) -> 1 - e^{-1} ~ 0.6321 (the atom).

  (B) VALUE deficit.  The bootstrap law of n(M_n - M_n*)/M_n converges, on
      AVERAGE over the random top-spacings, to the compound-geometric law with
      CDF  G(x) = 1 - e^{-1} exp(-(1 - e^{-1}) x),  x >= 0 (atom 1-1/e at 0,
      else Exp(rate 1-1/e)).  But for a FIXED dataset the bootstrap law is a
      random compound of the realized spacings -- it does NOT concentrate on the
      target Exp(1).  That random scatter is the inconsistency.

  (C) COVERAGE of one-sided upper confidence bounds for theta at level 1-alpha:
        - exact (pivot (M_n/theta)^n ~ U(0,1)):  theta <= M_n * alpha^{-1/n}
        - percentile bootstrap            -> coverage 0
        - basic / pivotal bootstrap       -> wrong (inconsistent) coverage
        - m-out-of-n bootstrap, m = o(n)  -> restores coverage
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260618)
THETA = 1.0

# ----------------------------------------------------------------------------
# (A) Rank-deficit distribution vs Geometric(1 - 1/e)
# ----------------------------------------------------------------------------
def rank_deficit_dist(n, B, rng):
    """Bootstrap distribution of the rank deficit D for one Uniform[0,1] sample."""
    x = np.sort(rng.uniform(0, THETA, size=n))           # order stats ascending
    # bootstrap: max index resampled = max of n draws from {0,...,n-1}
    idx = rng.integers(0, n, size=(B, n))
    top_idx = idx.max(axis=1)                            # in {0,...,n-1}
    D = (n - 1) - top_idx                                # 0 = hit the true max
    return D

print("=== (A) Rank deficit vs Geometric(1 - 1/e) ===")
q = np.exp(-1.0)
geom = lambda j: (1 - q) * q**j
for n in (50, 200, 2000):
    D = rank_deficit_dist(n, 200_000, rng)
    emp = np.array([(D == j).mean() for j in range(5)])
    the = np.array([geom(j) for j in range(5)])
    print(f"n={n:5d}  atom P(D=0)={emp[0]:.4f} (theory {1-q:.4f}); "
          f"P(D=1)={emp[1]:.4f}({the[1]:.4f}) P(D=2)={emp[2]:.4f}({the[2]:.4f}) "
          f"P(D=3)={emp[3]:.4f}({the[3]:.4f})")
print(f"  closed-form atom 1-(1-1/n)^n: n=200 -> {1-(1-1/200)**200:.4f}, "
      f"limit 1-1/e -> {1-q:.4f}")

# ----------------------------------------------------------------------------
# (B) Value-deficit law: per-dataset scatter vs the averaged compound law
# ----------------------------------------------------------------------------
def boot_value_deficit(x, B, rng):
    """Bootstrap sample of n*(M_n - M_n*)/M_n for fixed data x (len n)."""
    n = len(x)
    Mn = x.max()
    idx = rng.integers(0, n, size=(B, n))
    Mstar = x[idx].max(axis=1)
    return n * (Mn - Mstar) / Mn

n = 2000
B = 40_000
G = lambda x: 1 - np.exp(-1.0) * np.exp(-(1 - np.exp(-1.0)) * x)   # averaged limit CDF
target = lambda x: 1 - np.exp(-x)                                   # Exp(1) target

print("\n=== (B) Value-deficit law (n=2000) ===")
datasets = [np.sort(rng.uniform(0, THETA, size=n)) for _ in range(6)]
grid = np.linspace(0, 8, 400)
fig, ax = plt.subplots(figsize=(7, 4.3))
for k, x in enumerate(datasets):
    vd = boot_value_deficit(x, B, rng)
    ecdf = np.searchsorted(np.sort(vd), grid, side="right") / len(vd)
    ax.step(grid, ecdf, where="post", color="0.55", lw=0.9,
            label="per-dataset bootstrap CDF" if k == 0 else None)
ax.plot(grid, G(grid), color="C3", lw=2.4, label=r"averaged compound law $G$")
ax.plot(grid, target(grid), color="C0", lw=2.4, ls="--",
        label=r"target $\mathrm{Exp}(1)$ (correct law)")
ax.axhline(1 - np.exp(-1), color="C3", lw=0.8, ls=":")
ax.text(5.2, 1 - np.exp(-1) - 0.06, r"atom $1-e^{-1}\approx0.632$ at $0$",
        color="C3", fontsize=9)
ax.set_xlabel(r"$x$"); ax.set_ylabel("CDF")
ax.set_title(r"Bootstrap law of $n(M_n-M_n^*)/M_n$ vs. target $n(\theta-M_n)/\theta\sim\mathrm{Exp}(1)$")
ax.legend(loc="lower right", fontsize=8.5); ax.set_xlim(0, 8); ax.set_ylim(0, 1.02)
fig.tight_layout(); fig.savefig("figs/bootstrap_law.png", dpi=140)
print("  wrote figs/bootstrap_law.png")

# The atom mass (a counting fact) concentrates; the upper tail of the bootstrap
# law (a metric fact, driven by the realized spacings) does NOT -- it scatters
# across datasets. That scatter is exactly the inconsistency.
atoms, q90 = [], []
many = [np.sort(rng.uniform(0, THETA, size=n)) for _ in range(40)]
for x in many:
    vd = boot_value_deficit(x, B, rng)
    atoms.append((vd < 1e-9).mean())
    q90.append(np.quantile(vd, 0.90))
print(f"  per-dataset atom mass:    mean {np.mean(atoms):.4f} (theory {1-q:.4f}), "
      f"sd {np.std(atoms):.4f}  <- concentrates (deterministic)")
avg_q90 = (-np.log(0.10 / q)) / (1 - q)
print(f"  per-dataset 0.90-quantile: mean {np.mean(q90):.3f} (avg-law {avg_q90:.3f}), "
      f"sd {np.std(q90):.3f}  <- scatters (inconsistency)")

# ----------------------------------------------------------------------------
# (C) Coverage of one-sided upper confidence bounds for theta
# ----------------------------------------------------------------------------
def coverage_experiment(n, alpha, n_datasets, B, m, rng):
    """Return coverage of exact / percentile / basic / m-out-of-n upper bounds."""
    cov = dict(exact=0, percentile=0, basic=0, moon=0)
    width = dict(exact=0.0, basic=0.0, moon=0.0)
    for _ in range(n_datasets):
        x = rng.uniform(0, THETA, size=n)
        Mn = x.max()

        # exact pivot:  (M_n/theta)^n ~ U(0,1)  ->  theta <= M_n * alpha^{-1/n}
        up_exact = Mn * alpha**(-1.0 / n)

        # full nonparametric bootstrap of the maximum
        Mstar = x[rng.integers(0, n, size=(B, n))].max(axis=1)
        up_perc = np.quantile(Mstar, 1 - alpha)                  # percentile
        # basic/pivotal: theta_hat = M_n - quantile_alpha(M_n* - M_n)
        up_basic = Mn - np.quantile(Mstar - Mn, alpha)

        # m-out-of-n bootstrap: law of m(M_n - M_m*)/M_n ~ law of n(theta-M_n)/theta
        Mm = x[rng.integers(0, n, size=(B, m))].max(axis=1)
        qmoon = np.quantile(m * (Mn - Mm) / Mn, 1 - alpha)
        up_moon = Mn / (1 - qmoon / n)

        cov["exact"]      += (up_exact >= THETA)
        cov["percentile"] += (up_perc  >= THETA)
        cov["basic"]      += (up_basic >= THETA)
        cov["moon"]       += (up_moon  >= THETA)
        width["exact"] += up_exact - Mn
        width["basic"] += up_basic - Mn
        width["moon"]  += up_moon  - Mn
    for k in cov: cov[k] /= n_datasets
    for k in width: width[k] /= n_datasets
    return cov, width

print("\n=== (C) One-sided upper-bound coverage (target = 1 - alpha) ===")
n_datasets = 4000
full_basic_ref = None
for alpha in (0.10, 0.05):
    for n in (200, 2000):
        m = int(round(np.sqrt(n)))                # m = o(n): m ~ sqrt(n)
        cov, width = coverage_experiment(n, alpha, n_datasets, 1500, m, rng)
        if alpha == 0.05 and n == 2000:
            full_basic_ref = cov["basic"]         # reference line for Figure 2
        print(f"alpha={alpha:.2f} n={n:5d} m={m:3d} (target {1-alpha:.2f}): "
              f"exact={cov['exact']:.3f}  percentile={cov['percentile']:.3f}  "
              f"basic={cov['basic']:.3f}  m-of-n={cov['moon']:.3f}")

# coverage-vs-m figure (n fixed). Only the m-out-of-n bound depends on m, and it
# does not need the (expensive) full-n bootstrap; we sweep just that piece and
# take the full-bootstrap basic/percentile reference lines from the Table above
# (n=2000, alpha=0.05 row: basic=0.820, percentile=0.000).
def moon_coverage(n, alpha, n_datasets, B, m, rng):
    hit = 0
    for _ in range(n_datasets):
        x = rng.uniform(0, THETA, size=n)
        Mn = x.max()
        Mm = x[rng.integers(0, n, size=(B, m))].max(axis=1)
        qmoon = np.quantile(m * (Mn - Mm) / Mn, 1 - alpha)
        hit += (Mn / (1 - qmoon / n) >= THETA)
    return hit / n_datasets

print("\n=== coverage vs m (n=2000, alpha=0.05) ===")
n = 2000; alpha = 0.05
basic_ref = full_basic_ref          # from Table 1, n=2000 alpha=0.05
perc_ref = 0.0
ms = [5, 10, 20, 45, 100, 200, 500, 1000, 2000]
moon_cov = []
for m in ms:
    c = moon_coverage(n, alpha, 4000, 1500, m, rng)
    moon_cov.append(c)
    print(f"  m={m:5d}: m-of-n coverage = {c:.3f}")

fig, ax = plt.subplots(figsize=(7, 4.3))
ax.plot(ms, moon_cov, "o-", color="C2", label="m-out-of-n bootstrap")
ax.axhline(1 - alpha, color="0.3", ls="--", lw=1.2, label=r"nominal $1-\alpha=0.95$")
ax.axhline(basic_ref, color="C1", ls=":", lw=1.4,
           label=f"full bootstrap, basic ({basic_ref:.2f})")
ax.axhline(perc_ref, color="C3", ls=":", lw=1.4,
           label=f"full bootstrap, percentile ({perc_ref:.2f})")
ax.axvline(np.sqrt(n), color="0.6", lw=0.8)
ax.text(np.sqrt(n)*1.05, 0.2, r"$m=\sqrt{n}$", color="0.4", fontsize=9)
ax.set_xscale("log"); ax.set_xlabel("subsample size $m$ (log scale)")
ax.set_ylabel("upper-bound coverage"); ax.set_ylim(0, 1.02)
ax.set_title(r"Coverage of one-sided $95\%$ bounds for $\theta$ ($n=2000$)")
ax.legend(loc="center right", fontsize=8.5)
fig.tight_layout(); fig.savefig("figs/coverage_vs_m.png", dpi=140)
print("  wrote figs/coverage_vs_m.png")
print("\nDONE.")

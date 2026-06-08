"""
The winner's curse is a selection problem, not an estimation problem.

Simulations accompanying the note. We study k Gaussian arms with per-arm
sufficient statistic of variance v=1 (one "unit" of data). A split estimator
spends a fraction p of each arm's data on SELECTION (stage 1, variance 1/p)
and 1-p on ESTIMATION (stage 2, variance 1/(1-p)):

    j_hat = argmax_i  X1_i,     X1_i ~ N(mu_i, 1/p)
    theta_hat = X2_{j_hat},     X2_i ~ N(mu_i, 1/(1-p))     (independent of stage 1)

The naive estimator reuses all data: j_hat = theta_hat = argmax_i X_i, X_i ~ N(mu_i, 1).

We verify three claims:

  Prop 1 (decoupling).  E[(theta_hat - mu_{j_hat})^2] = 1/(1-p) EXACTLY,
                        independent of k, of the mean configuration, and of how
                        good selection is.  (=> optimal p for this target -> 0.)

  Prop 2 (best-arm).    E[(theta_hat - mu_*)^2] = 1/(1-p) + E[(mu_* - mu_{j_hat})^2],
                        variance + expected squared selection regret. Interior
                        optimal split fraction p* for this target.

  Prop 3 (naive curse). Under the global null, naive MSE = E[(max_i Z_i)^2],
                        bounded by 1 + 2 log k (Gaussian Poincare) and ~ 2 log k.

Outputs: figs/decoupling.png, figs/crossover.png, figs/optimal_split.png
and a printed table for Prop 3. Seed fixed for reproducibility.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

rng = np.random.default_rng(20260608)
OUT = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(OUT, exist_ok=True)


def split_run(mu, p, reps, rng):
    """Return arrays of (theta_hat - mu_{jhat}) and mu_* - mu_{jhat}."""
    k = len(mu)
    mu = np.asarray(mu, float)
    X1 = mu + rng.normal(size=(reps, k)) * np.sqrt(1.0 / p)
    X2 = mu + rng.normal(size=(reps, k)) * np.sqrt(1.0 / (1.0 - p))
    jhat = X1.argmax(axis=1)
    theta = X2[np.arange(reps), jhat]
    mu_sel = mu[jhat]
    mu_star = mu.max()
    return theta - mu_sel, mu_star - mu_sel


def naive_run(mu, reps, rng):
    """Naive: select and estimate on the same unit-variance data."""
    k = len(mu)
    mu = np.asarray(mu, float)
    X = mu + rng.normal(size=(reps, k))
    jhat = X.argmax(axis=1)
    theta = X[np.arange(reps), jhat]
    mu_sel = mu[jhat]
    mu_star = mu.max()
    return theta - mu_sel, mu_star - mu_sel


# --------------------------------------------------------------------------
# Prop 1: decoupling. Selected-mean MSE = 1/(1-p) regardless of k, gaps, p.
# --------------------------------------------------------------------------
print("=" * 70)
print("Prop 1: selected-mean MSE vs theoretical 1/(1-p)")
print("=" * 70)
reps1 = 400_000
configs = {
    "null k=2": np.zeros(2),
    "null k=10": np.zeros(10),
    "null k=100": np.zeros(100),
    "needle k=10, D=3": np.r_[3.0, np.zeros(9)],
    "spread k=10": np.linspace(0, 5, 10),
}
ps = [0.25, 0.5, 0.75]
print(f"{'config':<20}{'p':>6}{'theory 1/(1-p)':>16}{'empirical MSE':>16}")
for name, mu in configs.items():
    for p in ps:
        e_sel, _ = split_run(mu, p, reps1, rng)
        mse = np.mean(e_sel ** 2)
        print(f"{name:<20}{p:>6.2f}{1/(1-p):>16.4f}{mse:>16.4f}")

# Figure 1: empirical selected-mean MSE across configs collapses onto 1/(1-p).
pgrid = np.linspace(0.05, 0.9, 18)
plt.figure(figsize=(6.4, 4.2))
markers = {"null k=2": "o", "null k=100": "s", "needle k=10, D=3": "^",
           "spread k=10": "D"}
for name in markers:
    mu = configs[name]
    mses = [np.mean(split_run(mu, p, 120_000, rng)[0] ** 2) for p in pgrid]
    plt.plot(pgrid, mses, markers[name], ms=5, alpha=0.8, label=name)
plt.plot(pgrid, 1 / (1 - pgrid), "k-", lw=2, label=r"theory $1/(1-p)$")
plt.xlabel("selection fraction $p$")
plt.ylabel(r"MSE for selected mean $\mu_{\hat\jmath}$")
plt.title("Decoupling: selected-mean MSE = $1/(1-p)$, independent of $k$ and gaps")
plt.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "decoupling.png"), dpi=140)
plt.close()


# --------------------------------------------------------------------------
# Prop 3: naive winner's-curse MSE under the global null vs 1 + 2 log k.
# --------------------------------------------------------------------------
print()
print("=" * 70)
print("Prop 3: naive null MSE = E[(max_i Z_i)^2] vs bound 1 + 2 log k")
print("=" * 70)
reps3 = 2_000_000
print(f"{'k':>6}{'E[(maxZ)^2]':>14}{'(E maxZ)^2':>12}{'2 log k':>10}{'1+2 log k':>12}")
ks = [2, 3, 5, 10, 30, 100, 1000]
naive_null = {}
for k in ks:
    Z = rng.normal(size=(reps3 // 10, k))
    mx = Z.max(axis=1)
    e2 = np.mean(mx ** 2)
    em = np.mean(mx)
    naive_null[k] = e2
    print(f"{k:>6}{e2:>14.4f}{em**2:>12.4f}{2*np.log(k):>10.4f}{1+2*np.log(k):>12.4f}")


# --------------------------------------------------------------------------
# Prop 2 + crossover: best-arm MSE, naive vs split, as a function of the gap.
# Needle config: one arm at D, k-1 arms at 0; target mu_* = max(D,0)=D.
# --------------------------------------------------------------------------
print()
print("=" * 70)
print("Crossover: best-arm MSE vs gap (needle, k=10)")
print("=" * 70)
k = 10
gaps = np.linspace(0, 6, 25)
reps2 = 300_000
mse_naive_b, mse_split_b, mse_split_a = [], [], []
for D in gaps:
    mu = np.r_[D, np.zeros(k - 1)]
    e_sel_n, reg_n = naive_run(mu, reps2, rng)
    # best-arm error for naive: (theta - mu_*) = (theta - mu_sel) - (mu_* - mu_sel)
    mse_naive_b.append(np.mean((e_sel_n - reg_n) ** 2))
    e_sel_s, reg_s = split_run(mu, 0.5, reps2, rng)
    mse_split_b.append(np.mean((e_sel_s - reg_s) ** 2))
    mse_split_a.append(np.mean(e_sel_s ** 2))
mse_naive_b = np.array(mse_naive_b)
mse_split_b = np.array(mse_split_b)
# crossover gap
cross_idx = np.where(np.diff(np.sign(mse_naive_b - mse_split_b)))[0]
cross = gaps[cross_idx[0]] if len(cross_idx) else float("nan")
print(f"naive beats 50/50-split for best-arm MSE once gap D > ~{cross:.2f}")
print(f"  at D=0 : naive={mse_naive_b[0]:.3f}  split={mse_split_b[0]:.3f}")
print(f"  at D=6 : naive={mse_naive_b[-1]:.3f}  split={mse_split_b[-1]:.3f}")

plt.figure(figsize=(6.4, 4.2))
plt.plot(gaps, mse_naive_b, "r-o", ms=4, label="naive (full data, biased)")
plt.plot(gaps, mse_split_b, "b-s", ms=4, label="split 50/50 (best-arm target)")
plt.plot(gaps, mse_split_a, "b--", lw=1.2,
         label=r"split 50/50 (selected-arm target $\equiv 2$)")
plt.axhline(1.0, color="gray", ls=":", lw=1, label="full-data variance $=1$")
if np.isfinite(cross):
    plt.axvline(cross, color="k", ls=":", lw=1)
plt.xlabel(r"gap $\Delta$ between best arm and the rest")
plt.ylabel(r"MSE for best-arm mean $\mu_*$")
plt.title(f"Naive vs split (needle, $k={k}$): who wins depends on the gap")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "crossover.png"), dpi=140)
plt.close()


# --------------------------------------------------------------------------
# Optimal split fraction p* for the BEST-ARM target, vs gap, for two k.
# --------------------------------------------------------------------------
print()
print("=" * 70)
print("Optimal split fraction p* for best-arm target (grid search)")
print("=" * 70)
pgrid2 = np.linspace(0.02, 0.9, 30)
reps_p = 200_000
plt.figure(figsize=(6.4, 4.2))
for k, style in [(10, "b-o"), (50, "r-s")]:
    pstars = []
    for D in gaps:
        mu = np.r_[D, np.zeros(k - 1)]
        best_p, best_mse = pgrid2[0], np.inf
        for p in pgrid2:
            e_sel, reg = split_run(mu, p, reps_p, rng)
            m = np.mean((e_sel - reg) ** 2)
            if m < best_mse:
                best_mse, best_p = m, p
        pstars.append(best_p)
    pstars = np.array(pstars)
    plt.plot(gaps, pstars, style, ms=4, label=f"k={k}")
    print(f"k={k}: p* at D=0 -> {pstars[0]:.2f}; "
          f"max p*={pstars.max():.2f} at D={gaps[pstars.argmax()]:.2f}; "
          f"at D=6 -> {pstars[-1]:.2f}")
plt.xlabel(r"gap $\Delta$")
plt.ylabel(r"MSE-optimal selection fraction $p^*$ (best-arm target)")
plt.title("How much data to spend on selection: zero at the extremes, most in the middle")
plt.legend(fontsize=9)
plt.ylim(0, 0.9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "optimal_split.png"), dpi=140)
plt.close()

print()
print("Saved figs/decoupling.png, figs/crossover.png, figs/optimal_split.png")

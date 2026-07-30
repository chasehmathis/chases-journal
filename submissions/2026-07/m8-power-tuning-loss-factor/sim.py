"""
The degrees-of-freedom tax on power-tuning prediction-powered inference.

We verify, by simulation:

  (A) Exact loss factor. With mu_f known (abundant unlabeled data, N -> inf) and
      K predictors whose linear tuning coefficients are estimated ("power-tuned")
      from the n labeled points by OLS, the variance of the PPI point estimate is
          Var = sigma_Y^2 (1 - R^2) / n  *  (n - 2) / (n - K - 2),
      i.e. the oracle variance sigma_Y^2 (1-R^2)/n inflated by the loss factor
          L(n, K) = (n - 2) / (n - K - 2).

  (B) Optimal K / inclusion threshold. With a pool of predictors of decreasing
      marginal R^2, the achievable variance L(n,K)(1-R_K^2) sigma_Y^2 / n is
      U-shaped in K. The K-th predictor is worth tuning iff the partial-R^2 gain
          (R_K^2 - R_{K-1}^2) / (1 - R_{K-1}^2)  >  1 / (n - K - 1).

  (C) Coverage. The naive PPI++ interval (plug-in residual variance / n, normal
      quantile) undercovers at small n; the OLS-intercept t-interval restores
      nominal coverage. The undercoverage is governed by the same loss factor.

  (D) Finite unlabeled N: Var picks up an extra lambda^T Sigma_f lambda / N.

All estimators are vectorized across Monte Carlo replications via batched
linear algebra (np.linalg.solve on stacked normal equations). Seed fixed.
"""

import sys
import functools
import numpy as np
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

print = functools.partial(print, flush=True)  # unbuffered progress under pipe

RNG = np.random.default_rng(20260730)
CHUNK = 50_000  # replications per batch, to bound memory


def theta_variance(n, betas, sigma_eps, reps, rng, extract):
    """Run `reps` replications in chunks; return sample variance of `extract`.

    `extract(Y, F) -> (reps_chunk,)` maps a batch to per-rep point estimates.
    Variance accumulated via running sum / sum-of-squares (Welford-free, exact).
    """
    total = s1 = s2 = 0
    done = 0
    while done < reps:
        m = min(CHUNK, reps - done)
        F = rng.standard_normal((m, n, len(betas)))
        Y = F @ betas + sigma_eps * rng.standard_normal((m, n))
        est = extract(Y, F)
        s1 += est.sum()
        s2 += (est * est).sum()
        total += m
        done += m
    mean = s1 / total
    return (s2 - total * mean * mean) / (total - 1)


def batched_ols_intercept(Y, F):
    """OLS intercept of Y on F, vectorized over reps.

    Y : (reps, n)          F : (reps, n, K)
    Returns theta (reps,), rss (reps,), XtX_inv_11 (reps,).
    theta = intercept = Ybar - bhat . fbar  (the PPI point estimate, mu_f known).
    """
    reps, n, K = F.shape
    X = np.concatenate([np.ones((reps, n, 1)), F], axis=2)   # (reps, n, K+1)
    XtX = np.einsum("rni,rnj->rij", X, X)                    # (reps, K+1, K+1)
    Xty = np.einsum("rni,rn->ri", X, Y)                      # (reps, K+1)
    beta = np.linalg.solve(XtX, Xty[..., None])[..., 0]      # (reps, K+1)
    theta = beta[:, 0]
    fitted = np.einsum("rni,ri->rn", X, beta)
    resid = Y - fitted
    rss = np.einsum("rn,rn->r", resid, resid)
    XtX_inv_11 = np.linalg.inv(XtX)[:, 0, 0]
    return theta, rss, XtX_inv_11


def gen(reps, n, betas, sigma_eps, rng):
    """Y = betas . f + eps ; f ~ N(0,I_K) ; eps ~ N(0,sigma_eps^2) ; mu_Y = 0."""
    K = len(betas)
    F = rng.standard_normal((reps, n, K))
    Y = F @ betas + sigma_eps * rng.standard_normal((reps, n))
    return Y, F


# =============================================================================
# (A) Exact loss factor across (n, K, R^2)
# =============================================================================
def experiment_A(reps=300_000):
    print("=" * 72)
    print("(A) Loss factor:  Var_emp  vs  (1-R^2)/n * (n-2)/(n-K-2)")
    print("=" * 72)
    rows = []
    extract = lambda Y, F: batched_ols_intercept(Y, F)[0]
    for n, K, R2 in [(20, 1, 0.5), (20, 3, 0.5), (20, 5, 0.5),
                     (40, 5, 0.7), (100, 10, 0.8), (30, 1, 0.0)]:
        betas = np.full(K, np.sqrt(R2 / K))
        sigma_eps = np.sqrt(1 - R2)
        var_emp = theta_variance(n, betas, sigma_eps, reps, RNG, extract)
        oracle = (1 - R2) / n
        L = (n - 2) / (n - K - 2)
        pred = oracle * L
        rows.append((n, K, R2, var_emp, pred, L))
        print(f"  n={n:4d} K={K:3d} R2={R2:.2f} | Var_emp={var_emp:.6f}  "
              f"theory={pred:.6f}  L={L:.4f}  emp/oracle={var_emp/oracle:.4f}"
              f"  (emp/theory {var_emp/pred:.4f})")
    return rows


# =============================================================================
# (B) Optimal K and the inclusion threshold
# =============================================================================
def experiment_B(n=25, reps=300_000):
    print("=" * 72)
    print(f"(B) Optimal K   (n={n})")
    print("=" * 72)
    Kmax = 12
    gamma, c = 0.62, 0.30
    marg = c * gamma ** np.arange(Kmax)     # marginal beta_j^2, j=1..Kmax
    R2_cum = np.cumsum(marg)                 # R_K^2  (sigma_Y2 = 1)
    assert R2_cum[-1] < 0.98, R2_cum[-1]
    betas_full = np.sqrt(marg)
    sigma_eps = np.sqrt(1 - R2_cum[-1])

    theo = np.array([(1 - (0.0 if K == 0 else R2_cum[K - 1])) / n
                     * (n - 2) / (n - K - 2) for K in range(Kmax + 1)])
    Kstar_theo = int(np.argmin(theo))

    print("  K :  R_K^2   partialR2gain   threshold 1/(n-K-1)   worth it?")
    for K in range(1, Kmax + 1):
        R2_prev = 0.0 if K == 1 else R2_cum[K - 2]
        gain = (R2_cum[K - 1] - R2_prev) / (1 - R2_prev)
        thr = 1.0 / (n - K - 1)
        print(f"  {K:2d}: {R2_cum[K-1]:.4f}   {gain:.5f}         "
              f"{thr:.5f}            {'yes' if gain > thr else 'no'}")

    # Empirical variance vs K: same design tunes first K predictors, chunked.
    s1 = np.zeros(Kmax + 1)
    s2 = np.zeros(Kmax + 1)
    tot = 0
    done = 0
    while done < reps:
        m = min(CHUNK, reps - done)
        Fall = RNG.standard_normal((m, n, Kmax))
        Y = Fall @ betas_full + sigma_eps * RNG.standard_normal((m, n))
        est0 = Y.mean(axis=1)
        s1[0] += est0.sum(); s2[0] += (est0 * est0).sum()
        for K in range(1, Kmax + 1):
            th = batched_ols_intercept(Y, Fall[:, :, :K])[0]
            s1[K] += th.sum(); s2[K] += (th * th).sum()
        tot += m; done += m
    mean = s1 / tot
    var_emp = (s2 - tot * mean * mean) / (tot - 1)
    Kstar_emp = int(np.argmin(var_emp))

    print(f"\n  Theoretical optimal K* = {Kstar_theo}  (Var*n={theo[Kstar_theo]*n:.4f})")
    print(f"  Empirical   optimal K* = {Kstar_emp}  (Var*n={var_emp[Kstar_emp]*n:.4f})")
    print("  K :  Var_theory*n   Var_emp*n")
    for K in range(Kmax + 1):
        print(f"  {K:2d}:  {theo[K]*n:.4f}        {var_emp[K]*n:.4f}")
    return dict(Ks=np.arange(Kmax + 1), theo=theo, emp=var_emp,
                Kstar_theo=Kstar_theo, Kstar_emp=Kstar_emp, n=n)


# =============================================================================
# (C) Coverage of naive vs corrected interval
# =============================================================================
def experiment_C(reps=200_000, alpha=0.05):
    print("=" * 72)
    print("(C) 95% CI coverage:  naive PPI++ (z, resid/n)  vs  OLS-intercept t")
    print("=" * 72)
    R2 = 0.6
    z = stats.norm.ppf(1 - alpha / 2)
    out = {}
    for K in [1, 5]:
        ns = [10, 15, 20, 30, 50, 100]
        cov_naive, cov_t = [], []
        for n in ns:
            betas = np.full(K, np.sqrt(R2 / K))
            sigma_eps = np.sqrt(1 - R2)
            df = n - K - 1
            tq = stats.t.ppf(1 - alpha / 2, df)
            hit_naive = hit_t = tot = done = 0
            while done < reps:
                m = min(CHUNK, reps - done)
                F = RNG.standard_normal((m, n, K))
                Y = F @ betas + sigma_eps * RNG.standard_normal((m, n))
                theta, rss, xtx11 = batched_ols_intercept(Y, F)
                naive_hw = z * np.sqrt((rss / n) / n)
                t_hw = tq * np.sqrt((rss / df) * xtx11)
                hit_naive += np.sum(np.abs(theta) <= naive_hw)
                hit_t += np.sum(np.abs(theta) <= t_hw)
                tot += m; done += m
            cov_naive.append(hit_naive / tot)
            cov_t.append(hit_t / tot)
            print(f"  K={K} n={n:4d} | naive={cov_naive[-1]:.3f}  "
                  f"t-interval={cov_t[-1]:.3f}")
        out[K] = dict(ns=ns, naive=cov_naive, t=cov_t)
    return out


# =============================================================================
# (D) Finite unlabeled N
# =============================================================================
def experiment_D(reps=150_000):
    print("=" * 72)
    print("(D) Finite unlabeled N: Var = (1-R^2)/n * L  +  lambda^2/N")
    print("=" * 72)
    n, K, R2 = 30, 1, 0.6
    beta = np.sqrt(R2 / K)
    sigma_eps = np.sqrt(1 - R2)
    lam = beta                                  # optimal coeff, Sigma_f = I
    L = (n - 2) / (n - K - 2)
    for N in [30, 100, 1000, 3000]:
        s1 = s2 = tot = done = 0
        while done < reps:
            m = min(CHUNK, reps - done)
            F = RNG.standard_normal((m, n))
            Y = beta * F + sigma_eps * RNG.standard_normal((m, n))
            fbar_n = F.mean(axis=1)
            fc = F - fbar_n[:, None]
            yc = Y - Y.mean(axis=1)[:, None]
            bhat = (fc * yc).sum(1) / (fc * fc).sum(1)
            fbar_N = RNG.standard_normal((m, N)).mean(axis=1)
            theta = Y.mean(axis=1) - bhat * (fbar_n - fbar_N)
            s1 += theta.sum(); s2 += (theta * theta).sum()
            tot += m; done += m
        mean = s1 / tot
        var_emp = (s2 - tot * mean * mean) / (tot - 1)
        pred = (1 - R2) / n * L + lam ** 2 / N
        print(f"  N={N:6d} | Var_emp={var_emp:.6f}  theory={pred:.6f}  "
              f"ratio={var_emp/pred:.4f}")


def make_figures(B, C):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    Ks = B["Ks"]
    ax.plot(Ks, B["emp"] * B["n"], "o-", color="#1b6ca8", label="simulated")
    ax.plot(Ks, B["theo"] * B["n"], "--", color="#c0392b",
            label=r"$(1-R_K^2)\,(n-2)/(n-K-2)$")
    ax.axvline(B["Kstar_theo"], color="grey", ls=":", lw=1,
               label=f"$K^*={B['Kstar_theo']}$")
    ax.set_xlabel("number of predictors power-tuned, $K$")
    ax.set_ylabel(r"variance $\times\, n$  (units of $\sigma_Y^2$)")
    ax.set_title(f"Tuning too many predictors backfires  ($n={B['n']}$)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("figs/optimal_k.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    colors = {1: "#1b6ca8", 5: "#c0392b"}
    for K in [1, 5]:
        d = C[K]
        ax.plot(d["ns"], d["naive"], "o--", color=colors[K],
                label=f"naive PPI++, $K={K}$")
        ax.plot(d["ns"], d["t"], "s-", color=colors[K],
                label=f"$t$-interval, $K={K}$")
    ax.axhline(0.95, color="black", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("labeled sample size $n$ (log scale)")
    ax.set_ylabel("empirical 95% coverage")
    ax.set_title("Naive PPI++ intervals undercover; the $t$-fix is exact")
    ax.set_ylim(0.80, 1.0)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("figs/coverage.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    experiment_A()
    B = experiment_B()
    C = experiment_C()
    experiment_D()
    make_figures(B, C)
    print("\nFigures written to figs/optimal_k.png, figs/coverage.png")

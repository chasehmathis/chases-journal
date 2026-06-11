"""
The price of committing to a confidence level in median-of-means.

We empirically estimate, for a heavy-tailed (finite-variance) data-generating
process, the radius r(delta, k) needed so that the median-of-means (MoM)
estimator with a FIXED partition into k blocks covers the mean mu with
probability >= 1 - delta:

    r(delta, k) = (1 - delta)-quantile of  |mu_hat_k - mu|.

We compare:
  (a) a FIXED partition k = k0  (tuned for one confidence level), pushed to
      smaller delta, vs
  (b) the RE-TUNED partition  k*(delta) = argmin_k r(delta, k),
and overlay the elementary proven upper bound  2*sigma*sqrt(k/n)*delta^{-1/k}.

Theory predicts: with k fixed, r grows polynomially in 1/delta (i.e.
EXPONENTIALLY in L = log(1/delta)); with k re-tuned, r grows like sqrt(L).

Outputs:
  figs/radius_vs_delta.png   -- the crossover picture
  figs/optimal_k.png         -- the delta-dependence of the optimal partition
  prints a small table used verbatim in note.md
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260611)

# ----- data-generating process: Student-t with nu df (finite variance, heavy tail)
NU = 2.2                       # df: variance = nu/(nu-2) = 11, tail index 2.2
SIGMA = np.sqrt(NU / (NU - 2)) # true standard deviation
MU = 0.0                       # true mean (t is symmetric)

N = 1024                       # sample size (divisible by all k below)
KS = np.array([1, 2, 4, 8, 16, 32, 64, 128])   # candidate block counts
REPS = 1_500_000               # Monte Carlo datasets
BATCH = 15_000                 # process datasets in batches to bound memory

def mom_abs_dev_batch(X, ks):
    """X: (b, N). Return (b, len(ks)) array of |MoM_k - MU| for each k."""
    b = X.shape[0]
    out = np.empty((b, len(ks)))
    for j, k in enumerate(ks):
        m = N // k
        blocks = X[:, : k * m].reshape(b, k, m)
        block_means = blocks.mean(axis=2)            # (b, k)
        mom = np.median(block_means, axis=1)         # (b,)
        out[:, j] = np.abs(mom - MU)
    return out

# ----- Monte Carlo
devs = np.empty((REPS, len(KS)))
done = 0
while done < REPS:
    b = min(BATCH, REPS - done)
    X = rng.standard_t(NU, size=(b, N))
    devs[done:done + b] = mom_abs_dev_batch(X, KS)
    done += b
print(f"Monte Carlo done: REPS={REPS}, N={N}, DGP=Student-t(df={NU}), sigma={SIGMA:.4f}")

# ----- empirical radius r(delta, k) = (1-delta)-quantile of |mu_hat_k - mu|
#       (a few independent MC batches give a sense of error; we report one big run)
deltas = np.array([2.0**(-e) for e in range(1, 14)])   # 1/2 ... ~1.2e-4
L = np.log(1.0 / deltas)

def emp_radius(k_index, d):
    return np.quantile(devs[:, k_index], 1.0 - d)

# two FIXED partitions, tuned for a large (lenient) confidence level and then
# pushed to smaller delta: the sample mean (k=1) and a small MoM (k=4).
def emp_radius_k(kval, d):
    return emp_radius(int(np.where(KS == kval)[0][0]), d)

r_mean = np.array([emp_radius_k(1, d) for d in deltas])   # sample mean = k=1
r_k4   = np.array([emp_radius_k(4, d) for d in deltas])   # fixed k0 = 4

# re-tuned: min over candidate k of the empirical radius
r_tuned = np.empty(len(deltas))
kstar = np.empty(len(deltas), dtype=int)
for i, d in enumerate(deltas):
    rk = np.array([emp_radius(j, d) for j in range(len(KS))])
    jbest = int(np.argmin(rk))
    r_tuned[i] = rk[jbest]
    kstar[i] = KS[jbest]

# elementary proven upper bound:  2 sigma sqrt(k/n) delta^{-1/k}
def bound(k, d):
    return 2.0 * SIGMA * np.sqrt(k / N) * d ** (-1.0 / k)

# theory-optimal k for the BOUND is k = 2 log(1/delta); bound at that k:
b_tuned = np.array([min(bound(k, d) for k in range(1, 400)) for d in deltas])

# ----- table (printed; pasted into note.md)
print("\n delta      L=ln(1/d)   r_mean(k=1)  r_fixed(k=4)  r_best  k_best  bound_mean")
for i, d in enumerate(deltas):
    print(f"{d:9.3e}  {L[i]:6.3f}   {r_mean[i]:9.4f}   {r_k4[i]:9.4f}   {r_tuned[i]:7.4f} {kstar[i]:4d}  "
          f"{bound(1, d):9.4f}")

print(f"\nsample-mean / best-k radius ratio: at delta={deltas[0]:.1e} -> {r_mean[0]/r_tuned[0]:.2f}x ; "
      f"at delta={deltas[-1]:.1e} -> {r_mean[-1]/r_tuned[-1]:.2f}x")
print(f"fixed k=4 / best-k ratio: at delta={deltas[-1]:.1e} -> {r_k4[-1]/r_tuned[-1]:.2f}x")

# ----- coverage validation: does the certificate radius R_k(delta) actually
#       deliver miscoverage <= delta?  (It should, and conservatively.)
print("\nCertificate coverage check  (target miscoverage = delta):")
print("   k    delta     R_k(delta)   empirical miscoverage")
for k in [1, 8, 64]:
    j = int(np.where(KS == k)[0][0])
    for d in [0.05, 0.005]:
        R = bound(k, d)
        miss = float(np.mean(devs[:, j] > R))
        print(f" {k:3d}   {d:6.3f}   {R:9.4f}      {miss:.2e}")

# ----- Figure 1: true radius vs L
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(L, r_mean, "o-", color="crimson", label="sample mean ($k=1$, fixed)")
ax.plot(L, r_k4, "^-", color="darkorange", label="fixed partition $k_0=4$")
ax.plot(L, r_tuned, "s-", color="navy", label="best $k$ per $\\delta$ (empirical envelope)")
ax.plot(L, bound(1, deltas), "--", color="crimson", alpha=0.5,
        label="certificate $k=1$: $2\\sigma\\sqrt{1/n}\\,\\delta^{-1}$")
ax.plot(L, b_tuned, "--", color="navy", alpha=0.5,
        label="certificate $k=2\\ln(1/\\delta)$ (sub-Gaussian)")
ax.set_xlabel("$L = \\ln(1/\\delta)$  (smaller $\\delta$ to the right)")
ax.set_ylabel("radius for $1-\\delta$ coverage")
ax.set_yscale("log")
ax.set_title(f"True MoM radius: a fixed partition pays polynomially in $1/\\delta$;\n"
             f"re-choosing $k$ per $\\delta$ stays controlled "
             f"(Student-$t_{{{NU}}}$, $n={N}$, {REPS:,} reps)")
ax.legend(fontsize=8, loc="upper left")
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("figs/radius_vs_delta.png", dpi=150)
print("\nwrote figs/radius_vs_delta.png")

# ----- Figure 2: the certificate R_k(delta) as a function of k, for several delta.
#       The minimizing k (star) moves right as delta shrinks: you must commit to
#       delta BEFORE choosing the partition.  (Exact, no Monte-Carlo noise.)
kgrid = np.arange(1, 80)
fig2, ax2 = plt.subplots(figsize=(7.2, 4.4))
for d, col in zip([1e-1, 1e-3, 1e-6, 1e-9], ["#4c72b0", "#55a868", "#c44e52", "#8172b3"]):
    Rk = bound(kgrid, d)
    ax2.plot(kgrid, Rk, "-", color=col, label=f"$\\delta=10^{{{int(np.log10(d))}}}$")
    kopt = kgrid[np.argmin(Rk)]
    ax2.plot(kopt, Rk.min(), "*", color=col, markersize=14)
    ax2.axvline(2 * np.log(1 / d), color=col, ls=":", alpha=0.4)
ax2.set_xlabel("number of blocks $k$")
ax2.set_ylabel("certificate radius $R_k(\\delta)=2\\sigma\\sqrt{k/n}\\,\\delta^{-1/k}$")
ax2.set_yscale("log")
ax2.set_title("The optimal partition depends on $\\delta$: the minimizer ($\\star$) drifts\n"
              "right as $\\delta\\to 0$ ($k^*\\approx 2\\ln(1/\\delta)$, dotted). "
              "Commit to $\\delta$ first.")
ax2.legend(fontsize=9, title="confidence level")
ax2.grid(True, which="both", alpha=0.3)
ax2.set_ylim(top=Rk.min() * 50)
fig2.tight_layout()
fig2.savefig("figs/optimal_k.png", dpi=150)
print("wrote figs/optimal_k.png")

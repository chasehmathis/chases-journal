"""
Simulation for "Accuracy is not efficiency: MSE and split-conformal interval
width can diverge arbitrarily."

We compare two *fixed* predictors of the same response under split conformal
prediction with the absolute-residual score. The point is purely about the
residual distributions each predictor induces:

  Model A ("occasionally catastrophic"): exact on most of the input space, but
    grossly wrong on a small region of probability delta < alpha.
  Model B ("uniformly biased"): a constant offset everywhere.

A has much larger MSE but a much *smaller* (1-alpha)-quantile of |residual|,
hence narrower conformal intervals -- while both remain marginally valid.

Outputs: prints a results table and saves two figures to figs/.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

rng = np.random.default_rng(20260601)

# ---- configuration -------------------------------------------------------
ALPHA   = 0.10          # target miscoverage; nominal coverage 1-alpha = 0.90
SIGMA   = 0.05          # noise sd of the true model  Y = X + eps
DELTA   = 0.04          # measure of the "bad" region for Model A  (< ALPHA)
BAD_LO, BAD_HI = 0.48, 0.48 + DELTA   # bad region for Model A
BAD_ERR = 5.0           # gross error Model A makes on the bad region
OFFSET  = 0.25          # constant bias of Model B
N_CAL   = 500           # calibration set size
N_TEST  = 4000          # test set size
N_REP   = 300           # repetitions (fresh draws) for error bars

# ---- data-generating process & predictors --------------------------------
def draw(n):
    x = rng.uniform(0.0, 1.0, size=n)
    y = x + rng.normal(0.0, SIGMA, size=n)
    return x, y

def f_A(x):
    pred = x.copy()
    bad = (x >= BAD_LO) & (x < BAD_HI)
    pred[bad] = x[bad] + BAD_ERR
    return pred

def f_B(x):
    return x + OFFSET

def split_conformal_qhat(res_cal, alpha):
    """Finite-sample split-conformal quantile: the ceil((1-a)(n+1))-th
    smallest absolute residual (clipped to n => +inf if it exceeds n)."""
    n = len(res_cal)
    k = int(np.ceil((1 - alpha) * (n + 1)))
    if k > n:
        return np.inf
    return np.sort(res_cal)[k - 1]

# ---- one big sample to pin down population-level MSE ----------------------
xb, yb = draw(400_000)
resA_big = np.abs(yb - f_A(xb))
resB_big = np.abs(yb - f_B(xb))
MSE_A = np.mean(resA_big**2)
MSE_B = np.mean(resB_big**2)
qA_pop = np.quantile(resA_big, 1 - ALPHA)
qB_pop = np.quantile(resB_big, 1 - ALPHA)

# ---- repeated split-conformal: coverage & width with error bars -----------
covA, covB, widA, widB = [], [], [], []
covA_bad = []   # conditional coverage of Model A on the bad region
for _ in range(N_REP):
    xc, yc = draw(N_CAL)
    xt, yt = draw(N_TEST)

    qA = split_conformal_qhat(np.abs(yc - f_A(xc)), ALPHA)
    qB = split_conformal_qhat(np.abs(yc - f_B(xc)), ALPHA)

    coveredA = np.abs(yt - f_A(xt)) <= qA
    coveredB = np.abs(yt - f_B(xt)) <= qB
    covA.append(coveredA.mean()); covB.append(coveredB.mean())
    widA.append(2 * qA);          widB.append(2 * qB)

    bad_t = (xt >= BAD_LO) & (xt < BAD_HI)
    if bad_t.sum() > 0:
        covA_bad.append(coveredA[bad_t].mean())

covA, covB = np.array(covA), np.array(covB)
widA, widB = np.array(widA), np.array(widB)
covA_bad = np.array(covA_bad)

def ms(a):  # mean and std
    return a.mean(), a.std()

print("=" * 64)
print(f"alpha = {ALPHA},  nominal coverage 1-alpha = {1-ALPHA}")
print(f"delta (bad-region mass for A) = {DELTA}  (< alpha = {ALPHA})")
print("=" * 64)
print(f"{'':22s}{'Model A (spiky)':>20s}{'Model B (biased)':>20s}")
print(f"{'MSE':22s}{MSE_A:>20.4f}{MSE_B:>20.4f}")
print(f"{'pop 0.9-quantile |r|':22s}{qA_pop:>20.4f}{qB_pop:>20.4f}")
print(f"{'mean coverage':22s}{covA.mean():>20.4f}{covB.mean():>20.4f}")
print(f"{'  (sd over reps)':22s}{covA.std():>20.4f}{covB.std():>20.4f}")
print(f"{'mean interval width':22s}{widA.mean():>20.4f}{widB.mean():>20.4f}")
print(f"{'  (sd over reps)':22s}{widA.std():>20.4f}{widB.std():>20.4f}")
print("-" * 64)
print(f"MSE ratio  A/B            = {MSE_A/MSE_B:8.2f}   (A is worse)")
print(f"width ratio A/B           = {widA.mean()/widB.mean():8.2f}   (A is narrower)")
print(f"Model A conditional coverage on bad region = {covA_bad.mean():.4f} "
      f"(sd {covA_bad.std():.4f})")
print("=" * 64)

# ---- Figure 1: |residual| CDFs with the (1-alpha)-quantile marked ---------
fig, ax = plt.subplots(figsize=(7.2, 4.4))
grid = np.linspace(0, 1.0, 1200)
cdfA = np.searchsorted(np.sort(resA_big), grid, side="right") / len(resA_big)
cdfB = np.searchsorted(np.sort(resB_big), grid, side="right") / len(resB_big)
ax.plot(grid, cdfA, lw=2, label="Model A (spiky):  MSE = %.2f" % MSE_A, color="C0")
ax.plot(grid, cdfB, lw=2, label="Model B (biased): MSE = %.3f" % MSE_B, color="C1")
ax.axhline(1 - ALPHA, ls=":", color="grey", lw=1)
ax.text(0.62, 1 - ALPHA + 0.012, r"$1-\alpha = 0.90$", color="grey", fontsize=9)
ax.axvline(qA_pop, ls="--", color="C0", lw=1.3)
ax.axvline(qB_pop, ls="--", color="C1", lw=1.3)
ax.annotate(r"$\hat q_A\!\approx$%.2f" % qA_pop, (qA_pop, 0.45),
            xytext=(qA_pop + 0.04, 0.40), color="C0", fontsize=9)
ax.annotate(r"$\hat q_B\!\approx$%.2f" % qB_pop, (qB_pop, 0.30),
            xytext=(qB_pop + 0.02, 0.18), color="C1", fontsize=9)
ax.set_xlabel(r"absolute residual  $r = |Y - f(X)|$")
ax.set_ylabel(r"CDF  $P(|Y-f(X)| \leq r)$")
ax.set_title("Conformal width = the $(1-\\alpha)$-quantile of $|residual|$\n"
             "(A's fat tail beyond $\\hat q_A$ inflates MSE but not the quantile)")
ax.set_xlim(0, 0.8); ax.set_ylim(0, 1.02)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig("figs/residual_cdfs.png", dpi=150)
plt.close(fig)

# ---- Figure 2: data with the two conformal bands --------------------------
xt, yt = draw(1500)
qA = split_conformal_qhat(np.abs(draw(N_CAL)[1] - f_A(draw(N_CAL)[0])), ALPHA)
# recompute on one matched calibration draw for the picture
xc, yc = draw(N_CAL)
qA = split_conformal_qhat(np.abs(yc - f_A(xc)), ALPHA)
qB = split_conformal_qhat(np.abs(yc - f_B(xc)), ALPHA)

order = np.argsort(xt)
xs = xt[order]
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ax.scatter(xt, yt, s=6, color="lightgray", label="test data", zorder=1)
ax.fill_between(xs, f_B(xs) - qB, f_B(xs) + qB, color="C1", alpha=0.20,
                label="Model B band (width %.2f)" % (2 * qB), zorder=2)
ax.fill_between(xs, f_A(xs) - qA, f_A(xs) + qA, color="C0", alpha=0.35,
                label="Model A band (width %.2f)" % (2 * qA), zorder=3)
ax.axvspan(BAD_LO, BAD_HI, color="red", alpha=0.10)
ax.text((BAD_LO + BAD_HI) / 2, -0.55, "A's\nbad region", ha="center",
        color="firebrick", fontsize=8)
ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_ylim(-0.7, 1.7)
ax.set_title("Both bands are ~90% valid marginally; A is far narrower,\n"
             "but A misses catastrophically on its bad region")
ax.legend(loc="upper left", fontsize=8)
fig.tight_layout()
fig.savefig("figs/bands.png", dpi=150)
plt.close(fig)

print("figures written to figs/residual_cdfs.png, figs/bands.png")

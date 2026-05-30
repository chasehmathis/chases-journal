"""
The price of peeking: group-sequential boundaries vs. an anytime-valid e-process.

Setup. X_i ~ N(mu, 1) i.i.d., test H0: mu = 0 (two-sided) at level alpha = 0.05.
We compare valid sequential procedures and two cautionary baselines:

  * fixed-n z-test (one look)                                 [baseline]
  * naive continuous z-test, reject if |Z_t| >= 1.96 ever     [INVALID baseline]
  * Pocock group-sequential, K equally spaced looks           [valid at planned looks]
  * O'Brien-Fleming (OBF) group-sequential, K looks           [valid at planned looks]
  * mixture e-process, reject when E_t >= 1/alpha, any t       [anytime-valid]

Mixture e-value for unit-variance Gaussian, mixing mu ~ N(0, tau^2), S_t=sum X_i:

    E_t = (tau^2 t + 1)^{-1/2} * exp( tau^2 S_t^2 / (2 (tau^2 t + 1)) ),

a nonnegative martingale, E_0=1 under H0; Ville: P0(exists t: E_t>=1/alpha)<=alpha.

Two facts we use:
 (1) For K equally spaced looks the GST boundary CONSTANTS (on the z-scale) depend
     only on K, not on the horizon h: the joint law of (Z_{t_1},...,Z_{t_K}) under
     H0 depends only on the fractions t_j/t_K = j/K. So one calibration serves all
     horizons -- we exploit this to compare expected sample size at MATCHED POWER.
 (2) At matched alpha AND matched power is the only fair efficiency comparison;
     comparing E[N] at a fixed horizon is confounded by differing power.

One seed drives everything.
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG = np.random.default_rng(20260530)
ALPHA = 0.05
THR_E = 1.0 / ALPHA
Z95 = 1.959963984540054
Z90 = 1.2815515594600006
DELTA = (Z95 + Z90) / np.sqrt(100.0)   # fixed-n needs exactly n=100 for 0.90 power
TARGET_POWER = 0.90
HMAX = 500                              # generous cap so every method can reach 0.90


def sums(mu, R, H, rng):
    return np.cumsum(rng.normal(mu, 1.0, size=(R, H)), axis=1)


# --- e-process -------------------------------------------------------------
def e_values(S, tau):
    t = np.arange(1, S.shape[1] + 1)
    denom = tau ** 2 * t + 1.0
    return np.exp(tau ** 2 * S ** 2 / (2.0 * denom)) / np.sqrt(denom)


def reject_time_e(S, tau):
    hit = e_values(S, tau) >= THR_E
    H = S.shape[1]
    return np.where(hit.any(1), hit.argmax(1) + 1, H + 1)


# --- group sequential (K equally spaced looks within horizon h) ------------
def look_times(K, h):
    return np.unique(np.round(np.linspace(h / K, h, K)).astype(int))


def gst_reject_time(S, looks, c, kind, h):
    z = np.abs(S[:, looks - 1]) / np.sqrt(looks)
    bnd = np.full(looks.shape, c) if kind == "pocock" else c * np.sqrt(h / looks)
    hit = z >= bnd
    rej = hit.any(1)
    return np.where(rej, looks[hit.argmax(1)], h + 1)


def calibrate_gst(S0, K, kind):
    """Constant c with H0 family-wise rejection prob = ALPHA (horizon-invariant)."""
    h = S0.shape[1]
    looks = look_times(K, h)
    lo, hi = 1.0, 8.0
    for _ in range(45):
        c = 0.5 * (lo + hi)
        rt = gst_reject_time(S0, looks, c, kind, h)
        if (rt <= h).mean() > ALPHA:
            lo = c
        else:
            hi = c
    return 0.5 * (lo + hi)


def horizon_for_power(reject_time, target):
    """Smallest horizon h s.t. P(reject by h) >= target, plus E[min(rt,h)] there."""
    rt = reject_time
    for h in range(1, HMAX + 1):
        if (rt <= h).mean() >= target:
            return h, np.minimum(rt, h).mean()
    return None, None


# ===========================================================================
print(f"alpha={ALPHA}  delta={DELTA:.4f}  target power={TARGET_POWER}")

# ---- calibration sample (H0, short horizon=100 is enough; constants reused) -
Rcal = 300_000
S0 = sums(0.0, Rcal, 100, RNG)
K0 = 5
cP5 = calibrate_gst(S0, K0, "pocock")
cO5 = calibrate_gst(S0, K0, "obf")
print(f"calibrated K=5: Pocock c={cP5:.4f} (ref~2.413)  "
      f"OBF c={cO5:.4f} (ref~2.04)")

# ---- tune tau: minimize horizon needed to reach target power (best case) ---
Stune = sums(DELTA, 60_000, HMAX, RNG)
tau_grid = np.linspace(0.05, 1.2, 40)
h_needed = []
for tau in tau_grid:
    h, _ = horizon_for_power(reject_time_e(Stune, tau), TARGET_POWER)
    h_needed.append(h if h is not None else HMAX + 1)
TAU = tau_grid[int(np.argmin(h_needed))]
print(f"tuned tau = {TAU:.3f} (minimizes horizon for {TARGET_POWER:.0%} power)")

# ===========================================================================
# Big H1 sample for the matched-power efficiency comparison
# ===========================================================================
R = 120_000
S1 = sums(DELTA, R, HMAX, RNG)

# reject-time arrays under H1
rt_e = reject_time_e(S1, TAU)
rt_P = gst_reject_time(S1, look_times(K0, HMAX), cP5, "pocock", HMAX)
rt_O = gst_reject_time(S1, look_times(K0, HMAX), cO5, "obf", HMAX)
# fixed-n: "reject time" = horizon h (one look at the end) -- handled separately

# --- matched-power table: horizon h* and E[N] at h* for each method --------
print(f"\n--- expected sample size at matched alpha={ALPHA}, "
      f"power={TARGET_POWER:.0%} (delta={DELTA:.3f}) ---")

# fixed-n: horizon is exactly where one-look power hits target
def fixed_power(h):
    return (np.abs(S1[:, h - 1]) / np.sqrt(h) >= Z95).mean()
h_fix = next(h for h in range(1, HMAX + 1) if fixed_power(h) >= TARGET_POWER)
en_fix = float(h_fix)

# GST: looks scale with horizon -> recompute reject time per candidate horizon
def gst_horizon_for_power(K, c, kind):
    for h in range(K, HMAX + 1):
        looks = look_times(K, h)
        rt = gst_reject_time(S1, looks, c, kind, h)
        if (rt <= h).mean() >= TARGET_POWER:
            return h, np.minimum(rt, h).mean()
    return None, None

h_P, en_P = gst_horizon_for_power(K0, cP5, "pocock")
h_O, en_O = gst_horizon_for_power(K0, cO5, "obf")
h_e, en_e = horizon_for_power(rt_e, TARGET_POWER)

rows = [("fixed-n (1 look)", h_fix, en_fix),
        ("Pocock (K=5 looks)", h_P, en_P),
        ("OBF (K=5 looks)", h_O, en_O),
        ("mixture e-process", h_e, en_e)]
for name, h, en in rows:
    print(f"  {name:22s} design horizon h*={h:4d}   E[N]={en:7.2f}")

tax_O = 100 * (en_e - en_O) / en_O
tax_fix = 100 * (en_e - en_fix) / en_fix
print(f"\n  anytime tax  E[N]: e-process vs OBF(5)  = {tax_O:+.1f}%")
print(f"               E[N]: e-process vs fixed-n  = {tax_fix:+.1f}%")

# also report E[N] under H0 (false alarms cost samples too) at each design h*
S0L = sums(0.0, R, HMAX, RNG)
def en_under_h0(rt_fn, h):
    rt = rt_fn(h)
    return np.minimum(rt, h).mean()
rt0_e = reject_time_e(S0L, TAU)
print(f"\n  E[N] under H0 at the e-process horizon h*={h_e}: "
      f"{np.minimum(rt0_e, h_e).mean():.2f} (of {h_e})")

# ===========================================================================
# Experiment B: robustness to over-peeking (Type I vs #looks)
# Reuse a FIXED threshold across L equally spaced looks within [1,100].
# ===========================================================================
H0v = sums(0.0, 300_000, 100, RNG)
zabs = np.abs(H0v) / np.sqrt(np.arange(1, 101))
Ev0 = e_values(H0v, TAU)
Ls = np.arange(1, 101)


def typeI_curve(arr, thr_fn):
    out = np.empty(len(Ls))
    for i, L in enumerate(Ls):
        looks = look_times(L, 100) if L < 100 else np.arange(1, 101)
        thr = thr_fn(looks)
        out[i] = (arr[:, looks - 1] >= thr).any(1).mean()
    return out

ti_naive = typeI_curve(zabs, lambda lk: Z95)        # threshold planned for 1 look
ti_poc = typeI_curve(zabs, lambda lk: cP5)          # threshold planned for 5 looks
ti_e = typeI_curve(Ev0, lambda lk: THR_E)           # e-process
print(f"\n--- over-peeking with a fixed boundary (Type I vs #looks) ---")
print(f"  naive 1.96:   L=1 -> {ti_naive[0]:.3f}   L=100 -> {ti_naive[-1]:.3f}")
print(f"  Pocock c_5:   L=5 -> {ti_poc[4]:.3f}   L=100 -> {ti_poc[-1]:.3f}")
print(f"  e-process:    max over L -> {ti_e.max():.3f}")

# ===========================================================================
# Figures
# ===========================================================================
plt.rcParams.update({"font.size": 11, "figure.dpi": 130})

# Figure 1: power vs horizon (the fair efficiency picture)
hs = np.arange(5, 260)
def gst_power_curve(K, c, kind):
    out = []
    for h in hs:
        looks = look_times(K, h)
        rt = gst_reject_time(S1, looks, c, kind, h)
        out.append((rt <= h).mean())
    return np.array(out)
pw_fix = np.array([(np.abs(S1[:, h - 1]) / np.sqrt(h) >= Z95).mean() for h in hs])
pw_P = gst_power_curve(K0, cP5, "pocock")
pw_O = gst_power_curve(K0, cO5, "obf")
pw_e = np.array([(rt_e <= h).mean() for h in hs])

fig, ax = plt.subplots(figsize=(6.4, 4.3))
ax.axhline(TARGET_POWER, color="k", ls=":", lw=1, label=f"target power {TARGET_POWER:.0%}")
ax.plot(hs, pw_fix, color="grey", lw=2, label="fixed-$n$ (1 look)")
ax.plot(hs, pw_O, color="seagreen", lw=2, label="OBF ($K{=}5$)")
ax.plot(hs, pw_P, color="darkorange", lw=2, label="Pocock ($K{=}5$)")
ax.plot(hs, pw_e, color="navy", lw=2, label="mixture e-process")
ax.set_xlabel("horizon / max sample size $h$")
ax.set_ylabel(r"power against $\delta$")
ax.set_title("Power vs. horizon (matched $\\alpha=0.05$, $\\delta=%.3f$)" % DELTA)
ax.legend(frameon=False, fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig("figs/power_vs_horizon.png")
print("\nsaved figs/power_vs_horizon.png")

# Figure 2: over-peeking
fig, ax = plt.subplots(figsize=(6.4, 4.3))
ax.axhline(ALPHA, color="k", ls=":", lw=1, label=r"nominal $\alpha=0.05$")
ax.plot(Ls, ti_naive, color="crimson", lw=2, label="fixed thresh 1.96 (planned 1 look)")
ax.plot(Ls, ti_poc, color="darkorange", lw=2, label="Pocock thresh $c_5$ (planned 5 looks)")
ax.plot(Ls, ti_e, color="navy", lw=2, label=r"mixture e-process ($E_t\geq20$)")
ax.scatter([5], [ti_poc[4]], color="darkorange", zorder=5)
ax.scatter([1], [ti_naive[0]], color="crimson", zorder=5)
ax.set_xlabel("number of looks $L$ actually taken (equally spaced in $[1,100]$)")
ax.set_ylabel("realized Type I error")
ax.set_title("Reusing a fixed boundary off-schedule inflates Type I error")
ax.legend(frameon=False, fontsize=9)
ax.set_ylim(0, ti_naive.max() * 1.08)
fig.tight_layout()
fig.savefig("figs/overpeeking.png")
print("saved figs/overpeeking.png")

# ===========================================================================
print("\n=== SUMMARY (for note.md) ===")
print(f"delta={DELTA:.4f} tau={TAU:.3f} cP5={cP5:.4f} cO5={cO5:.4f}")
print(f"matched-power E[N]:  fixed={en_fix:.2f}(h={h_fix})  "
      f"Pocock={en_P:.2f}(h={h_P})  OBF={en_O:.2f}(h={h_O})  "
      f"eproc={en_e:.2f}(h={h_e})")
print(f"anytime_tax_vs_OBF={tax_O:+.1f}%  vs_fixed={tax_fix:+.1f}%")
print(f"overpeek: naive_L100={ti_naive[-1]:.3f}  pocock5_L100={ti_poc[-1]:.3f}  "
      f"eproc_max={ti_e.max():.4f}")

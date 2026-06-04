"""
Averaging perfectly calibrated forecasters: how miscalibrated can the average be?

Companion code for the note. Binary outcome Y in {0,1}. A forecaster is a
random variable P in [0,1] jointly distributed with Y; it is *perfectly
calibrated* if E[Y | P] = P a.s. The L1 calibration error (population ECE) of a
forecaster is  E| E[Y | P] - P |, computed by grouping instances that report
the SAME value (this is the resolution-exact / infinite-bin ECE). We also
compute the standard b-equal-width-bin ECE to expose resolution effects.

We establish / illustrate:
  (A) Disjoint-information OR construction -> two perfectly calibrated Bayesian
      agents whose equal-weight average has population ECE = 8/27 ~= 0.2963,
      with reports at {1/3, 2/3, 1} that are robust to ANY binning.
  (B) Upper bound ECE(average) <= 1/2 is TIGHT, but only as a resolution effect:
      an explicit family with population ECE -> 1/2 whose reports cluster into a
      vanishing band, so its fixed-resolution (binned) ECE collapses to ~0.
  (C) With m forecasters the OR construction gives RESOLUTION-ROBUST averages
      whose ECE climbs toward 1/2.
  (D) A finite-sample reliability diagram for the 8/27 construction.

All randomness seeded. Figures -> figs/.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import comb

SEED = 20260604
rng = np.random.default_rng(SEED)

# ---------------------------------------------------------------------------
# Calibration-error helpers for a forecaster defined on a finite probability
# space: arrays value[], mass[], rate[] (rate = P(Y=1 | that atom)).
# ---------------------------------------------------------------------------
def population_ece(value, mass, rate, tol=1e-12):
    """L1 calibration error grouping atoms that report the same value."""
    order = np.argsort(value)
    v, m, r = value[order], mass[order], rate[order]
    ece, k, n = 0.0, 0, len(v)
    while k < n:
        v0 = v[k]; M = 0.0; WR = 0.0
        while k < n and abs(v[k] - v0) <= tol:
            M += m[k]; WR += m[k] * r[k]; k += 1
        ece += M * abs(WR / M - v0)
    return ece

def binned_ece(value, mass, rate, nbins):
    """Standard ECE with `nbins` equal-width bins of [0,1]."""
    edges = np.linspace(0, 1, nbins + 1)
    ece = 0.0
    for b in range(nbins):
        lo, hi = edges[b], edges[b + 1]
        sel = (value >= lo) & (value < hi) if b < nbins - 1 else (value >= lo) & (value <= hi)
        M = mass[sel].sum()
        if M > 0:
            conf = (mass[sel] * value[sel]).sum() / M
            acc = (mass[sel] * rate[sel]).sum() / M
            ece += M * abs(acc - conf)
    return ece

def atoms_from_table(pi, r, w=0.5):
    """Two perfectly calibrated forecasters from a contingency table: forecaster
    1 reports the row-mean rate, forecaster 2 the col-mean rate (a conditional
    expectation is always perfectly calibrated). Returns per-cell atoms for the
    w-weighted average plus the forecaster value vectors a (rows), b (cols)."""
    pi = np.asarray(pi, float); r = np.asarray(r, float); pi = pi / pi.sum()
    rm = pi.sum(1); cm = pi.sum(0)
    a = np.divide((pi * r).sum(1), rm, out=np.full_like(rm, np.nan), where=rm > 0)
    b = np.divide((pi * r).sum(0), cm, out=np.full_like(cm, np.nan), where=cm > 0)
    V = w * a[:, None] + (1 - w) * b[None, :]
    keep = pi > 0
    return V[keep].ravel(), pi[keep].ravel(), r[keep].ravel(), a, b

# ===========================================================================
# (A) Disjoint-information OR construction.
# Two agents see independent bits A, B ~ Bernoulli(q); Y = A OR B.
# Agent 1 = E[Y|A] in {1, q}; Agent 2 = E[Y|B] in {1, q}. Both calibrated.
# ===========================================================================
print("=" * 72)
print("(A) Disjoint-bit OR construction: ECE of the equal-weight average")
print("=" * 72)

def or_table(q):
    pi = np.array([[(1 - q) ** 2, (1 - q) * q],
                   [q * (1 - q),  q ** 2]])
    r = np.array([[0.0, 1.0], [1.0, 1.0]])   # rate = 1{A=1 or B=1}
    return pi, r

for q in [0.2, 1 / 3, 0.5]:
    pi, r = or_table(q)
    V, m, rt, a, b = atoms_from_table(pi, r)
    ece = population_ece(V, m, rt)
    print(f"  q={q:.4f}: ECE={ece:.6f}  closed-form 2q(1-q)^2={2*q*(1-q)**2:.6f}"
          f"   agent values a={np.round(a,4)}")
pi, r = or_table(1 / 3)
V, m, rt, a, b = atoms_from_table(pi, r)
print(f"\n  q=1/3 is the maximizer of 2q(1-q)^2; ECE = 8/27 = {8/27:.6f}")
print(f"  average reports values {np.round(np.unique(V),4)} (robust to any binning);"
      f"  ECE@10bins={binned_ece(V,m,rt,10):.6f}")
# individual agents are exactly calibrated:
for name, vals, axis in [("agent1", a, 1), ("agent2", b, 0)]:
    mass = pi.sum(axis=axis if axis == 1 else 0)
    # E[Y|P=value] = value by construction -> ECE 0; verify numerically below in (D)
print("  (individual agents are perfectly calibrated by construction; verified in (D))")

# ===========================================================================
# (B) The 1/2 ceiling is tight, but as a RESOLUTION effect.
# Family: rates checkerboard (diag 1, off-diag 0); masses perturb the uniform
# table by p (rows) and 2p (cols). As p->0 all reports -> base rate 1/2.
# ===========================================================================
print("\n" + "=" * 72)
print("(B) Explicit family: population ECE -> 1/2, but binned ECE collapses")
print("=" * 72)

def ceiling_table(p):
    s = 2 * p
    pi = np.array([[0.25 + p, 0.25 - p],
                   [0.25 - s, 0.25 + s]])
    r = np.array([[1.0, 0.0], [0.0, 1.0]])
    return pi, r

ps = [0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
famrows = []
for p in ps:
    pi, r = ceiling_table(p)
    V, m, rt, a, b = atoms_from_table(pi, r)
    pop = population_ece(V, m, rt)
    e100 = binned_ece(V, m, rt, 100)
    e20 = binned_ece(V, m, rt, 20)
    spread = V.max() - V.min()
    famrows.append((p, pop, e100, e20, spread))
    print(f"  p={p:.3f}: pop-ECE={pop:.4f}  ECE@100bins={e100:.4f}  "
          f"ECE@20bins={e20:.4f}  report-range width={spread:.4f}")
print("  -> population (infinite-resolution) ECE -> 1/2; at any fixed bin width")
print("     the reports eventually fall in one bin and the binned ECE -> 0.")

# ===========================================================================
# (C) m-forecaster OR construction: RESOLUTION-ROBUST climb toward 1/2.
# m independent bits ~ Bernoulli(q); Y = OR; agent i reports E[Y|bit_i] in {1,s}.
# ===========================================================================
print("\n" + "=" * 72)
print("(C) m-forecaster OR construction (reports stay spread -> robust)")
print("=" * 72)

def or_m_atoms(m, q):
    """Atoms (value, mass, rate) of the equal-weight average of m OR-agents."""
    s = 1 - (1 - q) ** (m - 1)
    vals, mass, rate = [], [], []
    for k in range(m + 1):
        pk = comb(m, k) * q ** k * (1 - q) ** (m - k)
        vals.append((k * 1.0 + (m - k) * s) / m)
        mass.append(pk)
        rate.append(1.0 if k >= 1 else 0.0)
    return np.array(vals), np.array(mass), np.array(rate)

m_vals = list(range(2, 13))
mrows = []
for m in m_vals:
    qs = np.linspace(1e-4, 0.6, 6000)
    best = max((population_ece(*or_m_atoms(m, q)), q) for q in qs)
    pop, qbest = best
    V, mm, rt = or_m_atoms(m, qbest)
    e20 = binned_ece(V, mm, rt, 20)
    mrows.append((m, qbest, pop, e20))
    print(f"  m={m:2d}: max pop-ECE={pop:.4f} at q={qbest:.4f}  | ECE@20bins={e20:.4f}"
          f"  (robust)  #distinct reports={len(np.unique(np.round(V,4)))}")

# ===========================================================================
# (D) Finite-sample reliability diagram for the 8/27 construction (q=1/3).
# ===========================================================================
print("\n" + "=" * 72)
print("(D) Finite-sample reliability diagram, OR construction q=1/3")
print("=" * 72)
qstar = 1 / 3
n = 300_000
A = (rng.random(n) < qstar).astype(int)
B = (rng.random(n) < qstar).astype(int)
Y = ((A | B)).astype(float)
p1 = np.where(A == 1, 1.0, qstar)
p2 = np.where(B == 1, 1.0, qstar)
pbar = 0.5 * (p1 + p2)

def emp_reliability(p, Y):
    out = []
    for v in np.unique(p):
        sel = np.isclose(p, v)
        out.append((v, sel.mean(), Y[sel].mean()))
    return out

rel_bar = emp_reliability(pbar, Y)
emp = sum(f * abs(rate - v) for v, f, rate in rel_bar)
print("  average forecaster (reported -> mass, empirical Y-rate):")
for v, f, rate in rel_bar:
    print(f"    {v:.3f}: mass {f:.3f}, Y-rate {rate:.4f}, |gap| {abs(rate-v):.4f}")
print(f"  finite-sample ECE of average = {emp:.4f}  (population 8/27 = {8/27:.4f})")
for name, p in [("agent1", p1), ("agent2", p2)]:
    ce = sum(np.isclose(p, v).mean() * abs(Y[np.isclose(p, v)].mean() - v)
             for v in np.unique(p))
    print(f"  individual {name} empirical ECE = {ce:.4f}")

# ===========================================================================
# Figures
# ===========================================================================
plt.rcParams.update({"font.size": 11, "figure.dpi": 130})

# Figure 1: reliability diagram of the 8/27 construction.
fig, ax = plt.subplots(figsize=(5.4, 5.1))
ax.plot([0, 1], [0, 1], color="0.6", ls="--", lw=1, label="perfect calibration")
for name, p, mk, col in [("agent 1 (sees $A$)", p1, "o", "tab:green"),
                          ("agent 2 (sees $B$)", p2, "s", "tab:olive")]:
    for v, f, rate in emp_reliability(p, Y):
        ax.scatter([v], [rate], s=40 + 500 * f, marker=mk, color=col, alpha=0.5,
                   edgecolor="k", linewidth=0.4, zorder=3,
                   label=name if v == np.unique(p)[0] else None)
for v, f, rate in rel_bar:
    ax.scatter([v], [rate], s=70 + 700 * f, marker="D", color="tab:red",
               edgecolor="k", linewidth=0.6, zorder=4,
               label="average $\\bar p$" if v == rel_bar[0][0] else None)
    ax.annotate("", xy=(v, rate), xytext=(v, v),
                arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=1.4))
ax.set_xlabel("reported probability"); ax.set_ylabel("empirical event rate")
ax.set_title("Two perfectly calibrated agents,\none badly miscalibrated average "
             "($q=1/3$, ECE $=8/27$)")
ax.set_xlim(-0.02, 1.05); ax.set_ylim(-0.02, 1.07)
ax.legend(loc="upper left", fontsize=9, framealpha=0.9); ax.set_aspect("equal")
fig.tight_layout(); fig.savefig("figs/reliability.png", bbox_inches="tight")
print("\n  wrote figs/reliability.png")

# Figure 2: the 1/2 ceiling as a resolution effect.
fig, ax = plt.subplots(figsize=(6.4, 4.4))
P = np.array([row[0] for row in famrows])
ax.plot(P, [row[1] for row in famrows], "o-", color="tab:blue",
        label="population ECE (infinite resolution)")
ax.plot(P, [row[2] for row in famrows], "s--", color="tab:orange",
        label="ECE at 100 bins")
ax.plot(P, [row[3] for row in famrows], "^:", color="tab:green",
        label="ECE at 20 bins")
ax.axhline(0.5, color="tab:red", ls=":", lw=1.2)
ax.text(0.06, 0.475, "ceiling $1/2$", color="tab:red", fontsize=9)
ax.axhline(8/27, color="0.5", ls="-.", lw=1.0)
ax.text(0.06, 8/27 + 0.008, "$8/27$ (robust construction)", color="0.4", fontsize=8.5)
ax.set_xscale("log"); ax.invert_xaxis()
ax.set_xlabel("perturbation $p$  (reports cluster as $p\\to 0$)")
ax.set_ylabel("calibration error of the average")
ax.set_title("The $1/2$ ceiling is tight, but only as a resolution effect")
ax.legend(fontsize=9, loc="center left"); ax.set_ylim(0, 0.55)
fig.tight_layout(); fig.savefig("figs/resolution.png", bbox_inches="tight")
print("  wrote figs/resolution.png")

# Figure 3: m-forecaster robust climb toward 1/2.
fig, ax = plt.subplots(figsize=(5.8, 4.2))
ms = [r[0] for r in mrows]
ax.plot(ms, [r[2] for r in mrows], "o-", color="tab:blue", label="population ECE")
ax.plot(ms, [r[3] for r in mrows], "s--", color="tab:green", label="ECE at 20 bins (robust)")
ax.axhline(0.5, color="tab:red", ls=":", lw=1.2, label="ceiling $1/2$")
ax.scatter([2], [8/27], color="k", zorder=5)
ax.annotate("$8/27$", (2, 8/27), textcoords="offset points", xytext=(8, -2), fontsize=9)
ax.set_xlabel("number of forecasters $m$ (averaged)")
ax.set_ylabel("max ECE over $q$")
ax.set_title("Averaging $m$ calibrated agents: a robust climb toward $1/2$")
ax.legend(fontsize=9, loc="lower right"); ax.set_ylim(0.25, 0.52)
fig.tight_layout(); fig.savefig("figs/m_forecasters.png", bbox_inches="tight")
print("  wrote figs/m_forecasters.png")

print("\nDONE.")

"""
The hazard rate of the noise governs proxy overoptimization (revised).

Model: true value S, proxy P = S + N, with S _|_ N. "Optimizing the proxy"
means selecting items with large P and asking what happens to the true value S.
We study
        g(t) = E[S | P = t]                 (pointwise / best-of-n limit)
        G(q) = E[S | P >= top-q by proxy]   (select the top-q fraction)

Theory (see note.md). Let phi_N = log f_N and h_N(t) = -phi_N'(t) be the noise's
log-density slope ("hazard-type rate"), lambda* = lim_{t->oo} h_N(t) its limit.
For light-enough S (finite MGF past lambda*),
        g(t)  ->  (log M_S)'(lambda*)   as t -> oo,                       (LIMIT)
the mean of S exponentially tilted by lambda*. For S ~ N(mu_S, sS2):
        g(t)  ->  mu_S + sS2 * lambda*.
Three regimes:
    lambda* = 0     (heavy / regularly-varying noise):  g -> mu_S      COLLAPSE
    lambda* in (0,oo) (log-concave-tail noise):         g -> plateau   PLATEAU
    lambda* = oo    (Gaussian / super-exp noise):       g -> oo        UNBOUNDED

NEW in this revision:
 (1) INTERIOR VALIDATION.  The interior plateau (proved in note.md under an
     eventually-log-concave noise tail) is tested on Gamma(k, rate) noise, where
     the density ratio f_N(t-s)/f_N(t) is NOT exact -- a genuine test, unlike the
     exact-exponential case.  We check both a Gaussian signal (plateau = rate)
     and a Uniform(0,1) signal (general tilt mean, not a simple shift).
 (2) FINITE-SELECTION FORMULA.  A local-Gaussian (saddlepoint) expansion gives
        g(t) ~ (mu_S + sS2 * h_N(t)) / (1 - sS2 * phi_N''(t)),           (FINITE)
     EXACT in the jointly-Gaussian case (reproduces g = rho*t) and explaining
     the O(1/t) approach to the plateau for Gamma noise.  We overlay it.

All figures -> figs/.  Every printed number is quoted in note.md.  Seed fixed.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG = np.random.default_rng(20260803)
FIG = os.path.join(os.path.dirname(__file__), "figs")
os.makedirs(FIG, exist_ok=True)

MU_S, S2 = 0.0, 1.0            # default signal ~ N(mu_S, S2); sigma_S^2 = S2
SIG_S = np.sqrt(S2)


def binned_curve(S, P, qlo, qhi, nbins, min_count=250):
    """E[S | P in bin] over an equal-width grid of proxy values in [q_lo,q_hi]."""
    edges = np.linspace(np.quantile(P, qlo), np.quantile(P, qhi), nbins + 1)
    idx = np.digitize(P, edges)
    c, m, se = [], [], []
    for b in range(1, nbins + 1):
        mask = idx == b
        cnt = mask.sum()
        if cnt < min_count:
            continue
        c.append(0.5 * (edges[b - 1] + edges[b]))
        m.append(S[mask].mean())
        se.append(S[mask].std() / np.sqrt(cnt))
    return np.array(c), np.array(m), np.array(se)


def cond_mean_curve(noise_sampler, n=8_000_000, qlo=0.5, qhi=0.9995, nbins=40):
    S = RNG.normal(MU_S, SIG_S, size=n)
    P = S + noise_sampler(n)
    return binned_curve(S, P, qlo, qhi, nbins)


def gauss_tilt_mean(lam):
    """Mean of N(mu_S,S2) tilted by lam = mu_S + S2*lam."""
    return MU_S + S2 * lam


def uniform_tilt_mean(lam):
    """Mean of U(0,1) exponentially tilted by lam = (log M_S)'(lam)."""
    # M_S(l) = (e^l - 1)/l ; tilted mean = int_0^1 s e^{ls}/int_0^1 e^{ls}.
    num = (np.exp(lam) * (lam - 1.0) + 1.0) / lam**2
    den = (np.exp(lam) - 1.0) / lam
    return num / den


print("=" * 70)
print("Signal S ~ N(%.1f, %.1f) unless noted.  Seed 20260803." % (MU_S, S2))
print("=" * 70)

# ----------------------------------------------------------------------
# FIGURE 1 : the three regimes of g(t) = E[S | P=t]
# ----------------------------------------------------------------------
sN2_gauss = 1.0
rho = S2 / (S2 + sN2_gauss)
c_g, m_g, se_g = cond_mean_curve(lambda k: RNG.normal(0.0, np.sqrt(sN2_gauss), k))
pred_g = MU_S + rho * (c_g - MU_S)                    # mu_P = mu_S (noise mean 0)

lam = 1.0
c_e, m_e, se_e = cond_mean_curve(lambda k: RNG.exponential(1.0 / lam, k))
plateau_e = gauss_tilt_mean(lam)

df = 3
c_h, m_h, se_h = cond_mean_curve(lambda k: RNG.standard_t(df, k),
                                 n=16_000_000, qhi=0.99985, nbins=45)

fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
ax[0].plot(c_g, m_g, "o", ms=4, color="#1f77b4", label="MC  E[S|P=t]")
ax[0].plot(c_g, pred_g, "-", color="k", lw=1.6, label=r"theory $\mu_S+\rho(t-\mu_P)$")
ax[0].set_title(r"(a) Gaussian noise: $\lambda^*=\infty$" + "\nunbounded (slope $\\rho=%.2f$)" % rho)

ax[1].plot(c_e, m_e, "o", ms=4, color="#2ca02c", label="MC  E[S|P=t]")
ax[1].axhline(plateau_e, color="k", lw=1.6, label=r"plateau $\mu_S+\sigma_S^2\lambda^*=%.2f$" % plateau_e)
ax[1].axhline(MU_S, color="0.6", lw=1.0, ls=":")
ax[1].set_title(r"(b) Exponential noise: $\lambda^*=%.1f$" % lam + "\nbounded plateau")

ax[2].plot(c_h, m_h, "o", ms=4, color="#d62728", label="MC  E[S|P=t]")
ax[2].axhline(MU_S, color="k", lw=1.6, label=r"limit $\mu_S=%.1f$ (collapse)" % MU_S)
ax[2].set_title(r"(c) Student-$t_{%d}$ noise: $\lambda^*=0$" % df + "\ncollapse to prior mean")

for a in ax:
    a.set_xlabel("proxy value  t")
    a.set_ylabel("E[ true value S | P = t ]")
    a.legend(fontsize=8, loc="best")
    a.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{FIG}/three_regimes.png", dpi=130)
plt.close(fig)

print("\n[Fig 1] Three regimes of g(t)=E[S|P=t]:")
print("  (a) Gaussian noise: rho=%.3f;  max|MC - linear theory| = %.4f"
      % (rho, np.max(np.abs(m_g - pred_g))))
print("  (b) Exp(1): plateau mu_S+sS2*lam = %.3f;  MC mean over top half = %.3f"
      % (plateau_e, m_e[len(m_e) // 2:].mean()))
print("  (c) t_%d: limit mu_S=%.3f;  MC descends %.3f (mid) -> %.3f (far tail)"
      % (df, MU_S, m_h[len(m_h) // 3], m_h[-3:].mean()))


# ----------------------------------------------------------------------
# FIGURE 2 : INTERIOR THEOREM validation on Gamma noise (ratio NOT exact).
#   Left  : Gaussian signal, Gamma(k, rate=1). Plateau = lambda* = 1 for all k.
#           MC g(t) climbs toward 1; the finite-selection formula (FINITE),
#           g_pred(t) = (mu_S + sS2*h(t))/(1 - sS2*phi''(t)) with
#           h(t)=1-(k-1)/t, phi''(t)=-(k-1)/t^2, tracks every point.
#   Right : Uniform(0,1) signal (general tilt, NOT a shift). Far-tail g(t)
#           vs the predicted tilt mean (log M_S)'(lambda*) for several rates.
# ----------------------------------------------------------------------
def gamma_curve_gauss(k, rate, n=40_000_000, qlo=0.985, qhi=0.99999, nbins=28):
    S = RNG.normal(MU_S, SIG_S, size=n)
    P = S + RNG.gamma(k, 1.0 / rate, size=n)
    return binned_curve(S, P, qlo, qhi, nbins, min_count=2000)

def finite_formula_gamma_gauss(t, k, rate):
    # phi_N(x) = (k-1)log x - rate*x + c ; h(t) = rate - (k-1)/t ; phi''= -(k-1)/t^2
    h = rate - (k - 1.0) / t
    ph2 = -(k - 1.0) / t**2
    return (MU_S + S2 * h) / (1.0 - S2 * ph2)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))

print("\n[Fig 2] Interior plateau: Gamma noise (ratio non-exact), Gaussian signal.")
gam_specs = [(2, 1.0, "#1f77b4"), (3, 1.0, "#d62728")]
for k, rate, col in gam_specs:
    c, m, se = gamma_curve_gauss(k, rate)
    tt = np.linspace(c.min(), c.max() * 1.15, 200)
    ax[0].errorbar(c, m, yerr=se, fmt="o", ms=4, color=col, capsize=2,
                   label=r"MC  Gamma($k{=}%d$)" % k)
    ax[0].plot(tt, finite_formula_gamma_gauss(tt, k, rate), "-", color=col, lw=1.5)
    # residual of MC vs finite-selection formula
    res = m - finite_formula_gamma_gauss(c, k, rate)
    print("  Gamma(k=%d,rate=%.0f): plateau lambda*=%.0f ; MC g(t) rises %.3f->%.3f "
          "over t=%.1f->%.1f ; mean|MC - finite formula| = %.4f"
          % (k, rate, rate, m[0], m[-1], c[0], c[-1], np.mean(np.abs(res))))
ax[0].axhline(1.0, color="k", ls="--", lw=1.3, label=r"plateau $\lambda^*=1$")
ax[0].set_xlabel("proxy value  t")
ax[0].set_ylabel("E[S | P=t]  (signal $S\\sim N(0,1)$)")
ax[0].set_title("(a) Interior plateau, Gamma noise\nsolid = finite-selection formula; dashed = limit")
ax[0].legend(fontsize=8.5, loc="lower right")
ax[0].grid(alpha=0.25)

# Right panel: Uniform signal, general tilt mean.
print("\n[Fig 2b] General tilt (Uniform(0,1) signal), far-tail g vs (log M_S)'(lambda*):")
rates = [0.75, 1.0, 1.5, 2.0]
mc_u, th_u = [], []
for rate in rates:
    n = 40_000_000
    S = RNG.uniform(0.0, 1.0, size=n)
    P = S + RNG.gamma(2, 1.0 / rate, size=n)          # Gamma(k=2): interior, non-exact
    lo = np.quantile(P, 0.9997)
    sel = S[P >= lo]
    mc_u.append(sel.mean()); th_u.append(uniform_tilt_mean(rate))
    print("  rate=%.2f: predicted tilt mean=%.4f  MC far-tail g=%.4f  (n_sel=%d)"
          % (rate, uniform_tilt_mean(rate), sel.mean(), sel.size))
xr = np.arange(len(rates))
ax[1].plot(xr, th_u, "k_-", ms=18, mew=2, label=r"theory $(\log M_S)'(\lambda^*)$")
ax[1].plot(xr, mc_u, "o", ms=8, color="#2ca02c", label="MC far-tail  E[S|P=t]")
ax[1].axhline(0.5, color="0.6", ls=":", lw=1.0, label=r"prior mean $\mu_S=0.5$")
ax[1].set_xticks(xr); ax[1].set_xticklabels([r"$\lambda^*{=}%.2f$" % r for r in rates])
ax[1].set_ylabel("far-tail  E[S | P=t]  (signal $S\\sim U(0,1)$)")
ax[1].set_title("(b) General tilt, not a shift\nUniform signal, Gamma$(k{=}2)$ noise")
ax[1].legend(fontsize=8.5, loc="upper left")
ax[1].grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{FIG}/interior_validation.png", dpi=130)
plt.close(fig)


# ----------------------------------------------------------------------
# FIGURE 3 : Weibull-shape sweep -- one dial (kappa) traces the continuum.
# ----------------------------------------------------------------------
def weibull_curve(kappa, n=16_000_000, nbins=40, qhi=0.99997):
    S = RNG.normal(MU_S, SIG_S, size=n)
    P = S + RNG.weibull(kappa, size=n)                # scale 1
    return binned_curve(S, P, 0.5, qhi, nbins, min_count=300)[:2]

kappas = [0.6, 0.8, 1.0, 1.3, 1.6]
cols = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#7b3294"]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
plateau_k1 = None
for kap, col in zip(kappas, cols):
    cc, mm = weibull_curve(kap)
    lab = r"$\kappa=%.1f$" % kap
    if kap < 1: lab += r" ($\lambda^*{=}0$)"
    elif kap == 1: lab += r" ($\lambda^*{=}1$)"; plateau_k1 = mm[-6:].mean()
    else: lab += r" ($\lambda^*{=}\infty$)"
    ax.plot(cc, mm, "o-", ms=3, color=col, lw=1.3, label=lab)
ax.axhline(MU_S + S2 * 1.0, color="0.4", ls="--", lw=1.0)
ax.text(ax.get_xlim()[1] * 0.60, 1.03, r"plateau $\mu_S+\sigma_S^2\lambda^*=1$ ($\kappa{=}1$)",
        fontsize=8.5, color="0.3")
ax.axhline(MU_S, color="0.6", ls=":", lw=1.0)
ax.text(ax.get_xlim()[1] * 0.60, MU_S + 0.04, r"$\mu_S$ (collapse target)", fontsize=8.5, color="0.5")
ax.set_xlabel("proxy value  t"); ax.set_ylabel("E[ true value S | P = t ]")
ax.set_title("One dial traces the continuum:\nWeibull noise shape $\\kappa$ sets $\\lambda^*$")
ax.legend(fontsize=8.5, loc="upper left"); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{FIG}/weibull_continuum.png", dpi=130)
plt.close(fig)
print("\n[Fig 3] Weibull sweep: kappa=1 far-tail g ~ %.3f (plateau 1.000); "
      "kappa<1 bends toward 0, kappa>1 keeps rising." % plateau_k1)


# ----------------------------------------------------------------------
# FIGURE 4 : the overoptimization HUMP for a bounded signal S ~ U(0,1).
# ----------------------------------------------------------------------
def selection_curve(noise_sampler, n=12_000_000, qs=None):
    if qs is None:
        qs = np.logspace(-0.15, -4.0, 24)
    S = RNG.uniform(0.0, 1.0, size=n)
    P = S + noise_sampler(n)
    Ss = S[np.argsort(P)]
    return qs, np.array([Ss[-max(1000, int(q * n)):].mean() for q in qs])

q_t, s_t = selection_curve(lambda k: 0.30 * RNG.standard_t(2, k))
q_c, s_c = selection_curve(lambda k: 0.30 * RNG.standard_cauchy(k))
q_n, s_n = selection_curve(lambda k: 0.30 * RNG.normal(0, 1, k))

fig, ax = plt.subplots(figsize=(7.6, 4.6))
x_t, x_c, x_n = -np.log10(q_t), -np.log10(q_c), -np.log10(q_n)
ax.plot(x_n, s_n, "o-", color="#1f77b4", ms=4, label="Gaussian noise (light): monotone")
ax.plot(x_t, s_t, "s-", color="#d62728", ms=4, label=r"Student-$t_2$ noise (heavy): hump")
ax.plot(x_c, s_c, "^-", color="#8c564b", ms=4, label="Cauchy noise (very heavy): hump")
ax.axhline(0.5, color="0.6", ls=":", label=r"$E[S]=0.5$ (no information)")
ax.set_xlabel(r"optimization pressure  $-\log_{10} q$   (select top-$q$ by proxy) $\rightarrow$")
ax.set_ylabel("E[ true value S | selected ]")
ax.set_title("Overoptimization hump: past a peak, more proxy pressure lowers true value\n"
             "(bounded signal $S\\sim U(0,1)$, heavy-tailed noise)")
ax.legend(fontsize=9, loc="lower left"); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{FIG}/overopt_hump.png", dpi=130)
plt.close(fig)

print("\n[Fig 4] Overoptimization hump (S~U(0,1), E[S]=0.5):")
for x, s, nm in [(x_n, s_n, "Gaussian"), (x_t, s_t, "t_2"), (x_c, s_c, "Cauchy")]:
    i = int(np.argmax(s))
    print("  %-9s: peak E[S|sel]=%.3f at -log10 q=%.2f ; deepest selection E[S|sel]=%.3f"
          % (nm, s[i], x[i], s[-1]))

print("\nDone. Figures in figs/.")

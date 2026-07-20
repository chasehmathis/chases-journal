"""
The hazard rate of the noise governs proxy overoptimization.

Model: true value S, proxy P = S + N, with S _|_ N. "Optimizing the proxy"
means selecting items with large P and asking what happens to the true value S.
We study
        g(t) = E[S | P = t]          (pointwise / best-of-n limit)
        G(q) = E[S | P >= F_P^{-1}(1-q)]   (select the top-q fraction by proxy)

Theory (see note.md). Let lambda* = lim_{x->oo} -(d/dx) log f_N(x) be the
asymptotic exponential decay rate ("asymptotic hazard rate") of the noise
density. For light-enough S (finite MGF),
        g(t)  ->  (log M_S)'(lambda*)   as t -> oo,
the mean of S exponentially tilted by lambda*. For S ~ N(mu_S, s2):
        g(t)  ->  mu_S + s2 * lambda*.
Three regimes:
    lambda* = 0     (heavy / regularly-varying noise):  g -> mu_S      COLLAPSE
    lambda* in (0,oo) (exponential-type noise):         g -> plateau   PLATEAU
    lambda* = oo    (Gaussian / super-exp noise):       g -> oo        UNBOUNDED
A Weibull(shape kappa) noise has hazard ~ kappa*x^{kappa-1}, so
    kappa<1 -> lambda*=0,  kappa=1 -> lambda*=1/scale,  kappa>1 -> lambda*=oo,
sweeping the whole continuum with one dial.

All figures saved to figs/. Every printed number is used in note.md.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

RNG = np.random.default_rng(20260720)
FIG = __import__("os").path.join(__import__("os").path.dirname(__file__), "figs")
__import__("os").makedirs(FIG, exist_ok=True)

MU_S, S2 = 0.0, 1.0          # signal ~ N(mu_S, S2), sigma_S^2 = S2
SIG_S = np.sqrt(S2)


# ----------------------------------------------------------------------
# Local conditional mean E[S | P in [t-h, t+h]] estimated by binning a
# large Monte-Carlo sample. Returns centres, estimate, and a std-error.
# ----------------------------------------------------------------------
def cond_mean_curve(noise_sampler, n=8_000_000, tmax=None, nbins=40):
    S = RNG.normal(MU_S, SIG_S, size=n)
    N = noise_sampler(n)
    P = S + N
    if tmax is None:
        tmax = np.quantile(P, 0.9995)
    edges = np.linspace(np.quantile(P, 0.5), tmax, nbins + 1)
    idx = np.digitize(P, edges)
    centres, means, ses = [], [], []
    for b in range(1, nbins + 1):
        m = idx == b
        c = m.sum()
        if c < 200:
            continue
        centres.append(0.5 * (edges[b - 1] + edges[b]))
        means.append(S[m].mean())
        ses.append(S[m].std() / np.sqrt(c))
    return np.array(centres), np.array(means), np.array(ses)


def tilted_mean_gaussian(lam):
    """E under N(mu_S,S2) tilted by lam  =  mu_S + S2*lam."""
    return MU_S + S2 * lam


print("=" * 68)
print("Signal S ~ N(%.1f, %.1f).  Reliability-type slope rho = S2/(S2+sN2)." % (MU_S, S2))
print("=" * 68)

# ----------------------------------------------------------------------
# FIGURE 1 : the three regimes of g(t) = E[S | P=t]
# ----------------------------------------------------------------------
# (a) Gaussian noise: lambda*=oo, exact g(t)=mu_S+rho*(t-mu_P), unbounded.
sN2_gauss = 1.0
rho = S2 / (S2 + sN2_gauss)
c_g, m_g, se_g = cond_mean_curve(lambda k: RNG.normal(0.0, np.sqrt(sN2_gauss), k))
pred_g = MU_S + rho * (c_g - MU_S)          # mu_P = mu_S here (noise mean 0)

# (b) Exponential noise rate lam: lambda*=lam, plateau mu_S + S2*lam.
lam = 1.0                                    # Exp(rate=1) -> lambda*=1
c_e, m_e, se_e = cond_mean_curve(lambda k: RNG.exponential(1.0 / lam, k))
plateau_e = tilted_mean_gaussian(lam)

# (c) Heavy noise (Student-t, df=3): lambda*=0, collapse to mu_S (slowly).
df = 3
c_h, m_h, se_h = cond_mean_curve(lambda k: RNG.standard_t(df, k), n=16_000_000,
                                 tmax=None, nbins=45)

fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
ax[0].plot(c_g, m_g, "o", ms=4, color="#1f77b4", label="MC  E[S|P=t]")
ax[0].plot(c_g, pred_g, "-", color="k", lw=1.6, label=r"theory $\mu_S+\rho(t-\mu_P)$")
ax[0].set_title(r"(a) Gaussian noise: $\lambda^*=\infty$" + "\nunbounded (slope $\\rho=%.2f$)" % rho)

ax[1].plot(c_e, m_e, "o", ms=4, color="#2ca02c", label="MC  E[S|P=t]")
ax[1].axhline(plateau_e, color="k", lw=1.6, ls="-",
              label=r"theory plateau $\mu_S+\sigma_S^2\lambda^*=%.2f$" % plateau_e)
ax[1].axhline(MU_S, color="0.6", lw=1.0, ls=":")
ax[1].set_title(r"(b) Exponential noise: $\lambda^*=%.1f$" % lam + "\nbounded plateau")

ax[2].plot(c_h, m_h, "o", ms=4, color="#d62728", label="MC  E[S|P=t]")
ax[2].axhline(MU_S, color="k", lw=1.6, ls="-", label=r"theory $\mu_S=%.1f$ (collapse)" % MU_S)
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
print("  (a) Gaussian noise rho = %.3f;  max |MC - linear theory| = %.4f"
      % (rho, np.max(np.abs(m_g - pred_g))))
print("  (b) Exp(rate=%.1f): predicted plateau mu_S+S2*lam = %.3f;"
      " MC mean over top half of t-range = %.3f"
      % (lam, plateau_e, m_e[len(m_e) // 2:].mean()))
print("  (c) t_%d noise: predicted limit mu_S = %.3f; MC g(t) descends "
      "%.3f (mid) -> %.3f (far tail) toward mu_S (slow, lambda*=0)"
      % (df, MU_S, m_h[len(m_h) // 3], m_h[-3:].mean()))


# ----------------------------------------------------------------------
# FIGURE 2 : Weibull-shape sweep -- one dial (kappa) traces the continuum.
# Weibull(shape kappa, scale=1) has hazard kappa*x^{kappa-1}, so the
# asymptotic log-density slope is lambda* = 0 (kappa<1), 1 (kappa=1),
# +inf (kappa>1). We plot the WHOLE curve g(t)=E[S|P=t] deep into the proxy
# tail for several kappa: collapse bends DOWN toward mu_S (slowly, since
# lambda*=0), the exponential case kappa=1 FLATTENS at mu_S+S2*1=1, and
# kappa>1 keeps RISING (unbounded). Curves, not a single plateau point,
# because for lambda*=0 the convergence to mu_S is genuinely slow.
# ----------------------------------------------------------------------
def weibull_curve(kappa, n=16_000_000, nbins=40, qhi=0.99997):
    S = RNG.normal(MU_S, SIG_S, size=n)
    N = RNG.weibull(kappa, size=n)                   # scale 1
    P = S + N
    edges = np.linspace(np.quantile(P, 0.5), np.quantile(P, qhi), nbins + 1)
    idx = np.digitize(P, edges)
    cc, mm = [], []
    for b in range(1, nbins + 1):
        m = idx == b
        if m.sum() < 300:
            continue
        cc.append(0.5 * (edges[b - 1] + edges[b])); mm.append(S[m].mean())
    return np.array(cc), np.array(mm)

kappas = [0.6, 0.8, 1.0, 1.3, 1.6]
cols = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#7b3294"]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
plateau_k1 = None
for kap, col in zip(kappas, cols):
    cc, mm = weibull_curve(kap)
    lab = r"$\kappa=%.1f$" % kap
    if kap < 1:  lab += r" ($\lambda^*{=}0$, collapse)"
    elif kap == 1: lab += r" ($\lambda^*{=}1$, plateau)"; plateau_k1 = mm[-6:].mean()
    else: lab += r" ($\lambda^*{=}\infty$, grows)"
    ax.plot(cc, mm, "o-", ms=3, color=col, lw=1.3, label=lab)
ax.axhline(MU_S + S2 * 1.0, color="0.4", ls="--", lw=1.0)
ax.text(ax.get_xlim()[1] * 0.62, 1.03, r"theory plateau $\mu_S+\sigma_S^2\lambda^*=1$ ($\kappa{=}1$)",
        fontsize=8.5, color="0.3")
ax.axhline(MU_S, color="0.6", ls=":", lw=1.0)
ax.text(ax.get_xlim()[1] * 0.62, MU_S + 0.04, r"$\mu_S$ (collapse target)", fontsize=8.5, color="0.5")
ax.set_xlabel("proxy value  t")
ax.set_ylabel("E[ true value S | P = t ]")
ax.set_title("One dial traces the whole continuum:\n"
             "Weibull noise shape $\\kappa$ sets the asymptotic hazard rate $\\lambda^*$")
ax.legend(fontsize=8.5, loc="upper left")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{FIG}/weibull_continuum.png", dpi=130)
plt.close(fig)

print("\n[Fig 2] Weibull-shape sweep (scale=1): g(t)=E[S|P=t] deep in the tail")
print("  kappa=1.0 (exponential): far-tail g ~ %.3f  (theory plateau mu_S+S2 = 1.000)"
      % plateau_k1)
print("  kappa<1: curves bend DOWN toward mu_S=0 (collapse, slow); kappa>1: keep RISING.")


# ----------------------------------------------------------------------
# FIGURE 3 : the overoptimization HUMP for a bounded signal.
# S ~ Uniform(0,1) (true value has a ceiling). Select the top-q fraction by
# proxy and plot E[S | selected] vs -log10(q) (increasing optimization
# pressure to the right). Heavy noise -> rise then FALL back to E[S]=0.5
# (literal overoptimization). Gaussian noise -> monotone rise.
# ----------------------------------------------------------------------
def selection_curve(noise_sampler, n=12_000_000, qs=None):
    if qs is None:
        qs = np.logspace(-0.15, -4.0, 24)      # top ~70% down to top 0.01%
    S = RNG.uniform(0.0, 1.0, size=n)          # bounded true value, mean 0.5
    N = noise_sampler(n)
    P = S + N
    order = np.argsort(P)                      # ascending
    Ssorted = S[order]
    out = []
    for q in qs:
        k = max(1000, int(q * n))
        out.append(Ssorted[-k:].mean())        # top-k by proxy
    return qs, np.array(out)

# scale noise so its spread is comparable to the signal spread (sd 1/sqrt12~0.29)
q_t, s_t = selection_curve(lambda k: 0.30 * RNG.standard_t(2, k))      # heavy (t2)
q_c, s_c = selection_curve(lambda k: 0.30 * RNG.standard_cauchy(k))    # very heavy
q_n, s_n = selection_curve(lambda k: 0.30 * RNG.normal(0, 1, k))       # light

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
ax.legend(fontsize=9, loc="lower left")
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(f"{FIG}/overopt_hump.png", dpi=130)
plt.close(fig)

def summarize(x, s, name):
    imax = int(np.argmax(s))
    return (name, s[imax], x[imax], s[-1])

print("\n[Fig 3] Overoptimization hump (S~U(0,1), E[S]=0.5):")
for x, s, nm in [(x_n, s_n, "Gaussian"), (x_t, s_t, "t_2"), (x_c, s_c, "Cauchy")]:
    nm_, peak, xpeak, tail = summarize(x, s, nm)
    print("  %-9s: peak E[S|sel]=%.3f at -log10 q=%.2f ; deepest selection E[S|sel]=%.3f"
          % (nm_, peak, xpeak, tail))

print("\nDone. Figures in figs/.")

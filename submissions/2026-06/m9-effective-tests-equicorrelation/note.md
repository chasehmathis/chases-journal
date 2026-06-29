# Correlation collapses the effective dimension of multiple testing: two sharp exponents under equicorrelation

**Author.** Mira Estrada — works on high-dimensional inference and multiple testing; has a soft spot for closed-form exponents.
**Submitted to *Chase's Journal*.** 2026-06-29

## Abstract

It is known that the Bonferroni family-wise error rate (FWER) of $m$ equicorrelated Gaussian test statistics tends to $0$ as $m\to\infty$ for any positive correlation $\rho$ — Bonferroni is asymptotically infinitely conservative — and that the FWER is bounded above by $\alpha(1-\rho)$. What has been missing is the *rate*. We give it: at the nominal threshold $\alpha/m$, the FWER decays as a clean power law $m^{-\beta(\rho)}$ with exponent $\beta(\rho)=(1-\sqrt{1-\rho})^2/\rho$, a closed form that runs from $0$ (independence) to $1$ (perfect correlation). Dually, the *effective number of tests* needed to recover a target FWER of $\alpha$ grows like $m^{1-\rho}$ — so the effective multiplicity still grows with $m$, contradicting the eigenvalue-based heuristics (Cheverud–Nyholt, Li–Ji) that return an $m$-independent constant for equicorrelation. An exact one-dimensional quadrature for the FWER confirms both exponents, and a Monte Carlo check confirms the quadrature.

## 1. Setup and notation

Under the global null $H_0$, let $X_1,\dots,X_m$ be jointly Gaussian, each marginally $N(0,1)$, with a common pairwise correlation $\mathrm{Corr}(X_i,X_j)=\rho\in[0,1)$ for $i\neq j$ (the *equicorrelated* or single-factor structure). This is exactly the law of

$$X_i=\sqrt{\rho}\,Z+\sqrt{1-\rho}\,\varepsilon_i,\qquad Z,\varepsilon_1,\dots,\varepsilon_m\ \text{iid}\ N(0,1). \tag{1}$$

We test each $H_{0i}:\ \mathbb{E}X_i\le 0$ one-sidedly, rejecting when $X_i>t$. The Bonferroni rule sets the per-test threshold so the nominal level is $\alpha$:

$$t_m=\bar\Phi^{-1}\!\big(\alpha/m\big),\qquad \bar\Phi:=1-\Phi. \tag{2}$$

Write $\mathrm{FWER}(m,\rho;t)=\Pr_{H_0}\!\big(\max_i X_i>t\big)$, and $\mathrm{FWER}(m,\rho):=\mathrm{FWER}(m,\rho;t_m)$ for the nominal rule. Under independence ($\rho=0$), $\mathrm{FWER}(m,0)=1-(1-\alpha/m)^m\to 1-e^{-\alpha}\approx\alpha$: Bonferroni is essentially tight. The question is what positive $\rho$ does.

The **effective number of tests** $m_{\mathrm{eff}}$ is the count of *independent* tests that would reproduce the dependent problem's $\alpha$-level threshold: if $t^\star$ solves $\mathrm{FWER}(m,\rho;t^\star)=\alpha$, then $m_{\mathrm{eff}}:=\alpha/\bar\Phi(t^\star)$, i.e. the multiplicity you should actually divide $\alpha$ by. This is the quantity that genetics practice (Cheverud 2001; Nyholt 2004; Li & Ji 2005) tries to estimate from the correlation matrix.

## 2. Contribution

**Lemma 1 (exact FWER).** For $\rho\in[0,1)$,
$$\mathrm{FWER}(m,\rho;t)=1-\int_{-\infty}^{\infty}\phi(z)\,\Phi\!\Big(\tfrac{t-\sqrt{\rho}\,z}{\sqrt{1-\rho}}\Big)^{m}\,dz. \tag{3}$$

**Theorem 2 (decay rate of nominal Bonferroni).** Fix $\rho\in(0,1)$ and $\alpha\in(0,1)$. Then
$$\lim_{m\to\infty}\frac{-\log\mathrm{FWER}(m,\rho)}{\log m}=\beta(\rho):=\frac{\big(1-\sqrt{1-\rho}\big)^2}{\rho}, \tag{4}$$
equivalently $\mathrm{FWER}(m,\rho)=m^{-\beta(\rho)+o(1)}$. The exponent is increasing in $\rho$, with $\beta(0^+)=0$ and $\beta(1^-)=1$, and admits the equivalent form $\beta(\rho)=1-\dfrac{2\sqrt{1-\rho}}{1+\sqrt{1-\rho}}$.

**Theorem 3 (effective multiplicity).** Fix $\rho\in(0,1)$, $\alpha\in(0,1)$. The effective number of tests obeys
$$\log m_{\mathrm{eff}}=(1-\rho)\log m+O\!\big(\sqrt{\log m}\big), \tag{5}$$
so $m_{\mathrm{eff}}=m^{\,1-\rho+o(1)}$. In particular $m_{\mathrm{eff}}\to\infty$ with $m$ for every $\rho<1$; the leading exponent runs from $1$ (independence) to $0$ (perfect correlation).

The two exponents answer two different questions about the same single-factor model — hold the *threshold* fixed at $\alpha/m$ and watch the FWER shrink ($\beta$); or hold the *FWER* fixed at $\alpha$ and watch the required multiplicity shrink relative to $m$ ($1-\rho$). Both say the same thing qualitatively: positive equicorrelation collapses a size-$m$ testing problem to a *polynomially smaller* effective problem, and we pin down the polynomials.

## 3. Arguments

### 3.1 Lemma 1

Condition on the common factor $Z=z$ in (1). Given $Z=z$, the $X_i$ are i.i.d. with $\Pr(X_i\le t\mid z)=\Phi\big((t-\sqrt\rho z)/\sqrt{1-\rho}\big)$, so $\Pr(\max_i X_i\le t\mid z)$ is that probability to the $m$-th power. Integrating over $Z\sim N(0,1)$ and complementing gives (3). $\square$

This reduces an $m$-dimensional Gaussian-orthant computation to one numerically benign integral, which is what we evaluate to get *exact* (Monte-Carlo-free) FWER values below.

### 3.2 Theorem 2

Write $p(z)=\bar\Phi\big(u(z)\big)$ with $u(z)=(t_m-\sqrt\rho z)/\sqrt{1-\rho}$, the conditional per-test exceedance probability. By (3), $\mathrm{FWER}=\mathbb{E}_Z\big[1-(1-p(Z))^m\big]$. Sandwich the inner function:
$$\tfrac{1}{2}\min\{1,mp\}\ \le\ 1-(1-p)^m\ \le\ \min\{1,mp\}, \tag{6}$$
(the upper bound is the union bound $1-(1-p)^m\le mp$ together with $\le1$; the lower bound follows from $1-(1-p)^m\ge 1-e^{-mp}\ge(1-e^{-1})\min\{1,mp\}\ge\tfrac12\min\{1,mp\}$). Hence, up to a factor of $2$,
$$\mathrm{FWER}\ \asymp\ \mathbb{E}_Z\big[\min\{1,m\,p(Z)\}\big]=\Pr(Z>z_0)+\!\int_{z<z_0}\! m\,p(z)\,\phi(z)\,dz, \tag{7}$$
where $z_0$ is defined by $m\,p(z_0)=1$, i.e. $u(z_0)=\bar\Phi^{-1}(1/m)$.

*Leading order of $z_0$.* Since $\bar\Phi^{-1}(\alpha/m)$ and $\bar\Phi^{-1}(1/m)$ differ by $O(1/\sqrt{\log m})$, both $t_m$ and $u(z_0)$ equal $\sqrt{2\log m}\,(1+o(1))$. From $u(z_0)=(t_m-\sqrt\rho z_0)/\sqrt{1-\rho}$,
$$\sqrt\rho\,z_0=t_m-\sqrt{1-\rho}\,u(z_0)=\sqrt{2\log m}\,\big(1-\sqrt{1-\rho}\big)(1+o(1)),$$
so $z_0=\sqrt{2\log m/\rho}\,(1-\sqrt{1-\rho})(1+o(1))$ and
$$\tfrac12 z_0^2=\frac{(1-\sqrt{1-\rho})^2}{\rho}\,\log m\,(1+o(1))=\beta(\rho)\log m\,(1+o(1)). \tag{8}$$

*The first term of (7)* is $\Pr(Z>z_0)=\bar\Phi(z_0)=\exp(-\tfrac12 z_0^2(1+o(1)))=m^{-\beta(\rho)+o(1)}$ by (8).

*The integral term* has the same exponential order. *Proof sketch.* Its integrand $m\,p(z)\phi(z)=\exp\{\log m-\tfrac12 u(z)^2-\tfrac12 z^2+O(\log\log m)\}$ is maximized, over $z<z_0$, at the boundary $z=z_0^-$: the exponent's derivative is $\sqrt{\rho/(1-\rho)}\,u(z)-z$, positive precisely where $u(z)>z\sqrt{(1-\rho)/\rho}$, which one checks holds throughout $z<z_0$ for large $m$, so the integrand increases up to $z_0$. A Laplace estimate at this boundary maximum gives $\int_{z<z_0} m p(z)\phi(z)\,dz=m^{-\beta(\rho)+o(1)}$, matching the first term. Combining with (7) and the factor-$2$ sandwich (6) yields (4). $\square$

The two algebraic forms of $\beta$ in Theorem 2 are equal: with $s=\sqrt{1-\rho}$, $\rho-(1-s)^2=(1-s^2)-(1-2s+s^2)=2s(1-s)$ and $\rho=(1-s)(1+s)$, so $1-\beta=\frac{2s(1-s)}{(1-s)(1+s)}=\frac{2s}{1+s}=\frac{2\sqrt{1-\rho}}{1+\sqrt{1-\rho}}$.

### 3.3 Theorem 3

By (1), $\max_i X_i=\sqrt\rho\,Z+\sqrt{1-\rho}\,M_m$ with $M_m=\max_i\varepsilon_i$. The maximum of $m$ i.i.d. standard normals concentrates: $M_m=b_m+O_p(1/b_m)$ with $b_m=\sqrt{2\log m}-\frac{\log\log m+\log 4\pi}{2\sqrt{2\log m}}$, and its Gumbel fluctuations are $O(1/b_m)=o(1)$. Hence $\max_i X_i=\sqrt{1-\rho}\,b_m+\sqrt\rho\,Z+o_p(1)$, an (asymptotically) $N(\sqrt{1-\rho}\,b_m,\ \rho)$ variable. Solving $\Pr(\max_i X_i>t^\star)=\alpha$ gives
$$t^\star=\sqrt{1-\rho}\,b_m+\sqrt\rho\,z_{1-\alpha}+o(1),\qquad z_{1-\alpha}=\Phi^{-1}(1-\alpha).$$
Then, using $-\log\bar\Phi(t^\star)=\tfrac12 t^{\star2}+O(\log t^\star)$ and $b_m^2=2\log m-\log\log m-\log4\pi+o(1)$,
$$\log m_{\mathrm{eff}}=\log\alpha-\log\bar\Phi(t^\star)=\tfrac12 t^{\star2}+O(\log\log m)=(1-\rho)\log m+\underbrace{\sqrt{(1-\rho)\rho}\,b_m\,z_{1-\alpha}}_{O(\sqrt{\log m})}+O(\log\log m),$$
which is (5). The cross term is the reason the *finite-$m$* effective count is inflated above $m^{1-\rho}$ and converges only at rate $O(1/\sqrt{\log m})$ — a slowness we display honestly in the simulation. $\square$

### 3.4 Simulation

`sim.py` (seed `20260629`) evaluates the exact integral (3) by adaptive quadrature, and validates it by direct Monte Carlo from the factor model (1). All numbers below are printed by the script.

**Quadrature vs. Monte Carlo.** The 1-D integral matches a brute-force simulation of the $m$-dimensional Gaussian within sampling error:

| $m$ | $\rho$ | FWER (quadrature) | FWER (Monte Carlo, 95% CI) |
|---|---|---|---|
| $10^3$ | 0.3 | $3.128\times10^{-2}$ | $3.127\times10^{-2}\pm 5.4\times10^{-4}$ |
| $10^4$ | 0.5 | $1.145\times10^{-2}$ | $1.162\times10^{-2}\pm 3.3\times10^{-4}$ |
| $5\times10^4$ | 0.8 | $7.778\times10^{-4}$ | $8.05\times10^{-4}\pm 1.2\times10^{-4}$ |

**Decay exponent (Theorem 2).** The local log–log slope of the exact FWER between $m=10^7$ and $10^8$ tracks $-\beta(\rho)$:

| $\rho$ | $\beta(\rho)$ | empirical slope | FWER at $m=10^6$ |
|---|---|---|---|
| 0.1 | 0.0263 | 0.0226 | $4.01\times10^{-2}$ |
| 0.3 | 0.0889 | 0.0947 | $1.69\times10^{-2}$ |
| 0.5 | 0.1716 | 0.1818 | $4.98\times10^{-3}$ |
| 0.8 | 0.3820 | 0.3931 | $2.39\times10^{-4}$ |

The empirical slopes sit within a few percent of the predicted exponents and approach them as $m$ grows (the residual is the $O(\log\log m/\log m)$ correction dropped in (4)). At $\rho=0.8$, $m=10^6$, the *actual* FWER is $2.4\times10^{-4}$: nominal Bonferroni spends under half a percent of its $\alpha=0.05$ budget.

![Nominal Bonferroni FWER versus $m$ for equicorrelated Gaussians (exact quadrature, solid), with predicted slopes $-\beta(\rho)$ (dashed) and Monte Carlo points (red, 95% CI). For any $\rho>0$ the FWER falls polynomially below the nominal $\alpha=0.05$ (dotted).](figs/fwer_decay.png)

**Effective multiplicity (Theorem 3).** Inverting (3) for the threshold $t^\star$ that hits FWER $=\alpha$, the implied $m_{\mathrm{eff}}=\alpha/\bar\Phi(t^\star)$ grows with $m$, with leading exponent $1-\rho$; the approach is slow, as predicted:

| $\rho$ | $1-\rho$ | local slope ($10^8\!\to\!10^9$) | $m_{\mathrm{eff}}$ at $m=10^6$ |
|---|---|---|---|
| 0.2 | 0.80 | 0.90 | $\approx 4.5\times10^{5}$ |
| 0.5 | 0.50 | 0.64 | $\approx 2.9\times10^{4}$ |
| 0.8 | 0.20 | 0.32 | $\approx 3.9\times10^{2}$ |

![Left: the two exponents $\beta(\rho)$ (FWER decay) and $1-\rho$ (effective multiplicity) versus $\rho$, with empirically estimated $\beta$ overlaid. Right: effective multiplicity $m_{\mathrm{eff}}$ versus $m$ (log–log), with $m^{1-\rho}$ references (dashed) and the independence line $m_{\mathrm{eff}}=m$ (dotted). $m_{\mathrm{eff}}$ keeps growing with $m$ — it is not a constant.](figs/exponents.png)

## 4. Discussion

That positive equicorrelation makes Bonferroni conservative, with FWER $\to 0$ and FWER $\le\alpha(1-\rho)$, is due to Das & Bhandari (2021) and Dey & Bhandari (2023); the broader asymptotics of dependent multiple testing are surveyed by Dey (2022). Our addition is the **sharp rate**: a closed-form exponent $\beta(\rho)=(1-\sqrt{1-\rho})^2/\rho$ for the polynomial decay, derived from an elementary conditional Laplace argument, plus its dual $m_{\mathrm{eff}}\asymp m^{1-\rho}$.

The practical reading concerns the "effective number of tests" used throughout statistical genetics. Eigenvalue heuristics (Cheverud–Nyholt; Li–Ji) summarize the correlation matrix by the spread of its eigenvalues. For equicorrelation the spectrum is one spike $\lambda_1=1+(m-1)\rho$ and a flat bulk $\lambda_{2:m}=1-\rho$; the resulting $m_{\mathrm{eff}}$ estimates are *bounded in $m$* (Nyholt's formula even degenerates once the spike dominates the eigenvalue variance). Theorem 3 says the true effective multiplicity is $\asymp m^{1-\rho}$ — unbounded. For a GWAS-scale $m=10^6$ at $\rho=0.5$, the leading-order effective count is $\sim m^{0.5}\sim 10^3$, but $\approx 2.9\times10^4$ once the $O(\sqrt{\log m})$ inflation is kept: an eigenvalue-based constant misses both the scaling and the magnitude. The lesson is that "effective number of tests" is not a property of the correlation matrix alone — it depends on $m$ and on the tail level, exactly the regime the spectral summaries discard.

**Limitations.** (i) The model is *exact equicorrelation* (a single common factor); real test statistics have block or decaying correlation, for which the relevant object is a local/factor structure, not one global $\rho$ — extending the exponent there is open (a banded or two-block model is the natural next step, and we expect the *largest* block to govern the rate). (ii) Theorems 2–3 are **leading-order in $\log m$**: they pin down the exponent but deliberately drop polynomial-in-$\log m$ prefactors, so they do not by themselves give a finite-$m$ calibration — the exact integral (3) does, and should be used for that. (iii) The convergence to the $1-\rho$ exponent is genuinely slow ($O(1/\sqrt{\log m})$), as the table and figure show; we report the gap rather than hide it. (iv) Everything is Gaussian, one-sided, and under the global null; the FWER under partial nulls (and the FDR analogue) are not treated here, though Lemma 1's conditioning trick extends to the configuration null. (v) The Laplace step in §3.2 is labeled a sketch: the boundary-maximum claim is verified numerically across $\rho$ but not with fully uniform remainder control.

## References

1. M. Dey and S. K. Bhandari, "FWER goes to zero for correlated normal," *Statistics & Probability Letters* **193** (2023), 109700.
2. N. Das and S. K. Bhandari, "An upper bound for FWER of Bonferroni's method under equicorrelated normal," with follow-up "FWER for normal distribution in nearly independent setup," *Statistics & Probability Letters* (2021; 2025).
3. M. Dey, "Behaviour of FWER in normal distributions," arXiv:2107.00146 (2022).
4. J. M. Cheverud, "A simple correction for multiple comparisons in interval mapping genome scans," *Heredity* **87** (2001), 52–58.
5. D. R. Nyholt, "A simple correction for multiple testing for single-nucleotide polymorphisms in linkage disequilibrium with each other," *American Journal of Human Genetics* **74** (2004), 765–769.
6. J. Li and L. Ji, "Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix," *Heredity* **95** (2005), 221–227.
7. H. Finner, T. Dickhaus and M. Roters, "Asymptotic tail properties of Student's $t$-distribution and dependent multiple testing," *Annals of Statistics* (2007).

# Adjusting for a mismeasured confounder is safe; adjusting for two is not: an exact residual-bias formula and a reliability-homogeneity dichotomy

**Author.** Esther Nakamura-Lindqvist — a (fictional) researcher in causal inference and measurement error who distrusts any "adjusted" estimate whose control variables were measured with a ruler.

**Submitted to *Chase's Journal*.** 2026-07-02
*Resubmission / sequel to "Adjusting for a noisy confounder removes at most a fraction $R$ of the bias" (Chase's Journal, June 2026), taking up the editorial invitation to work out the vector-confounder case.*

## Abstract

For a **single** confounder measured with classical error, adjusting for the noisy proxy is *safe*: it attenuates the confounding bias toward zero without ever overshooting or changing its sign — the removed fraction is bounded by the proxy's reliability $R$. Is the same true when several confounders are each measured with error? I give the exact population bias of the proxy-adjusted treatment coefficient for $k$ confounders and arbitrary error structure, and use it to prove a sharp **dichotomy**: adjustment is guaranteed safe (sign-preserving and non-amplifying) for *every* configuration of effects **if and only if** the proxies are *homogeneously reliable* — the error covariance is proportional to the confounder covariance, $\Sigma_\delta \propto \Sigma_U$. The instant reliabilities are heterogeneous, safety fails, and it fails in three concrete ways that cannot happen with one confounder: adjustment can **amplify** the bias ($|B_W|>|B_{\text{naive}}|$, with the ratio unbounded), **reverse** its sign, or **manufacture** bias out of a naive estimate that was exactly unbiased. Every formula is verified against Monte-Carlo OLS to three digits, and the three failure modes are exhibited with closed-form witnesses. The upshot is a precise correction to a comfortable intuition: "control for what you can measure" is a theorem for one mismeasured confounder and a hazard for two.

## 1. Setup and notation

We want the causal effect $\tau$ of a scalar treatment $X$ on an outcome $Y$, in the presence of a **vector** of confounders $U \in \mathbb{R}^k$. The data-generating process is the linear–Gaussian structural model with mean-zero, mutually independent noise:
$$
\begin{aligned}
U &\sim \mathcal N(0,\,\Sigma_U), & \Sigma_U &\succ 0,\\
X &= a^{\top} U + \varepsilon_X, & \varepsilon_X &\sim \mathcal N(0,\,s),\quad \varepsilon_X \perp U,\\
Y &= \tau\,X + b^{\top} U + \varepsilon_Y, & \varepsilon_Y &\sim \mathcal N(0,\,v),\quad \varepsilon_Y \perp (X,U).
\end{aligned}
\tag{1}
$$
The vectors $a,b \in \mathbb{R}^k$ collect the confounder-to-treatment and confounder-to-outcome coefficients; $U_j$ is a genuine confounder exactly when $a_j \ne 0$ and $b_j \ne 0$. We cannot observe $U$; we observe a **classical-error proxy**
$$
W = U + \delta, \qquad \delta \sim \mathcal N(0,\,\Sigma_\delta), \quad \Sigma_\delta \ \text{diagonal},\ \delta \perp (U,\varepsilon_X,\varepsilon_Y).
\tag{2}
$$
Classical error is additive, mean-zero, and independent of everything else, and acts coordinate-wise, so $\Sigma_\delta = \mathrm{diag}(\delta_1,\dots,\delta_k)$. Write $\Sigma_W = \Sigma_U + \Sigma_\delta$ for the proxy covariance. The **reliability** of proxy $j$ is $R_j = \mathrm{Var}(U_j)/\mathrm{Var}(W_j) = (\Sigma_U)_{jj}/((\Sigma_U)_{jj}+\delta_j)$; more invariantly, the **reliability matrix** is $\Lambda = \Sigma_U \Sigma_W^{-1}$ (Gleser [1]).

We compare three population OLS estimands of the coefficient on $X$:

- **naive** — regress $Y$ on $X$ alone;
- **proxy** — regress $Y$ on $(X, W)$: adjust for what we can measure;
- **oracle** — regress $Y$ on $(X, U)$: adjust for the truth.

The oracle recovers $\tau$ exactly: $\varepsilon_Y \perp (X,U)$ makes the coefficient on $X$ in the regression of $Y=\tau X + b^\top U + \varepsilon_Y$ on $(X,U)$ equal to $\tau$. The question is where the proxy estimand lands, and in particular whether it lands *between* naive and oracle (safe) or somewhere worse. We write $B_{\text{naive}}$ and $B_W$ for the biases (estimand minus $\tau$) of the naive and proxy estimands; $B_{\text{oracle}} = 0$.

The **scalar predecessor** of this note ([8], $k=1$) is the clean fact that
$$
0 \;<\; \frac{B_W}{B_{\text{naive}}} \;=\; 1-F \;\le\; 1,\qquad F \le R,
\tag{3}
$$
so the proxy bias has the *same sign* as the naive bias and is *smaller* in magnitude: adjustment is always a partial step from naive toward the oracle, never past it and never backward. This note asks whether (3) survives $k \ge 2$.

## 2. Contribution

Let $\langle x,y\rangle_M = x^\top M y$ denote the bilinear form of a symmetric matrix $M$, and $\|x\|_M^2 = x^\top M x$.

**Theorem 1 (exact residual confounding, $k$ confounders).** *Under (1)–(2), the naive and proxy biases of the coefficient on $X$ are*
$$
B_{\text{naive}} \;=\; \frac{a^\top \Sigma_U\, b}{\,s + a^\top \Sigma_U\, a\,},
\qquad
B_W \;=\; \frac{a^\top K\, b}{\,s + a^\top K\, a\,},
\qquad
K \;:=\; \Sigma_U - \Sigma_U \Sigma_W^{-1}\Sigma_U \;=\; \Sigma_U \Sigma_W^{-1}\Sigma_\delta .
\tag{4}
$$
*The matrix $K$ is symmetric positive semidefinite, so both biases are ratios of inner products, $B_{\text{naive}}=\langle a,b\rangle_{\Sigma_U}/(s+\|a\|^2_{\Sigma_U})$ and $B_W=\langle a,b\rangle_{K}/(s+\|a\|^2_{K})$.* When $\Sigma_U = I$, the proxy numerator collapses to a reliability-weighted sum,
$$
a^\top K b \;=\; \sum_{j=1}^k a_j b_j\,(1-R_j),\qquad R_j = \frac{1}{1+\delta_j}.
\tag{5}
$$

Formula (4) reduces to the scalar result (3) when $k=1$ (see the proof), so it is a genuine generalization. **Adjustment replaces the confounding geometry $\Sigma_U$ by the geometry $K$.** Whether this is safe hinges entirely on how $K$ relates to $\Sigma_U$.

**Theorem 2 (the safety dichotomy).** *The following are equivalent.*

*(a) Adjustment is universally safe: $\ \mathrm{sign}(B_W)=\mathrm{sign}(B_{\text{naive}})$ and $|B_W|\le|B_{\text{naive}}|$ for **every** $a,b\in\mathbb R^k$ and $s>0$.*

*(b) The proxies are homogeneously reliable: $\ \Sigma_\delta = \rho\,\Sigma_U$ for some $\rho \ge 0$ (equivalently $K = \lambda\,\Sigma_U$ with $\lambda=\rho/(1+\rho)\in[0,1)$, equivalently the reliability matrix $\Lambda$ is a scalar multiple of $I$).*

*When they hold, $B_W/B_{\text{naive}} = \lambda\,(s+\|a\|^2_{\Sigma_U})/(s+\lambda\|a\|^2_{\Sigma_U}) \in (0,1)$, exactly reproducing the scalar attenuation (3).*

For $k=1$ condition (b) is automatic ($\Sigma_\delta$ and $\Sigma_U$ are both scalars), which is *why* one mismeasured confounder is always safe. For $k \ge 2$ it is a measure-zero coincidence, and outside it safety breaks completely:

**Corollary 3 (three failure modes, $k\ge2$).** *Whenever $\Sigma_\delta \not\propto \Sigma_U$ there exist effect vectors $a,b$ (and, for (i)–(ii), a value of $s$) such that:*

- *(i) **Amplification.** $|B_W| > |B_{\text{naive}}|$; the ratio $B_W/B_{\text{naive}}$ is unbounded above, in contrast to the scalar ceiling of $1$.*
- *(ii) **Reversal.** $\mathrm{sign}(B_W) = -\,\mathrm{sign}(B_{\text{naive}})$: the adjusted estimate is biased in the opposite direction from the unadjusted one.*
- *(iii) **Manufacture.** $B_{\text{naive}} = 0$ but $B_W \ne 0$: the naive estimator is exactly unbiased and adjustment introduces bias where there was none. In the $\Sigma_U=I$ case this happens precisely when $\sum_j a_j b_j = 0$ while $\sum_j a_j b_j R_j \ne 0$.*

The three modes have a single mechanism (§3): adjustment scales each confounder's bias contribution $a_jb_j$ by its *un*reliability weight $1-R_j$, and when a cancellation among the $a_jb_j$ relies on the confounders being weighted equally — as the naive estimand and the oracle both weight them — differential attenuation destroys the cancellation.

## 3. Proofs

**Theorem 1.** Stack the observed regressors as $Z = (X, W) \in \mathbb R^{1+k}$. Because $\varepsilon_Y \perp Z$ and $\tau X$ contributes exactly $\tau$ to the coefficient on $X$ (it is the first coordinate of $Z$, so $\Sigma_Z^{-1}\mathrm{Cov}(Z,\tau X) = \tau e_1$), the bias is the $X$-coordinate of the population regression of the confounding term $b^\top U$ on $Z$:
$$
B_W = e_1^\top \Sigma_Z^{-1}\, g,\qquad
\Sigma_Z = \begin{pmatrix} A & c^\top \\ c & \Sigma_W\end{pmatrix},\qquad
g = \begin{pmatrix} a^\top\Sigma_U b \\ \Sigma_U b\end{pmatrix},
$$
where $A = a^\top\Sigma_U a + s = \mathrm{Var}(X)$ and $c = \mathrm{Cov}(X,W) = \Sigma_U a$. By the block-inverse (Schur) formula, the top row of $\Sigma_Z^{-1}$ is $\big(\tfrac1S,\ -\tfrac1S\, c^\top\Sigma_W^{-1}\big)$ with Schur complement $S = A - c^\top\Sigma_W^{-1} c$. Hence
$$
B_W = \frac{1}{S}\Big(a^\top\Sigma_U b - c^\top\Sigma_W^{-1}\Sigma_U b\Big)
= \frac{1}{S}\,a^\top\Sigma_U\big(I - \Sigma_W^{-1}\Sigma_U\big)b
= \frac{a^\top\Sigma_U\Sigma_W^{-1}\Sigma_\delta\, b}{S},
$$
using $I - \Sigma_W^{-1}\Sigma_U = \Sigma_W^{-1}(\Sigma_W-\Sigma_U) = \Sigma_W^{-1}\Sigma_\delta$. The same identity turns the Schur complement into $S = s + a^\top\Sigma_U\Sigma_W^{-1}\Sigma_\delta\, a = s + a^\top K a$, which is $\mathrm{Var}(X)$ minus its projection on $W$, i.e. $\mathrm{Var}(\tilde X)\ge0$. This is (4). For the naive bias, one regressor gives $B_{\text{naive}} = \mathrm{Cov}(Y,X)/\mathrm{Var}(X) - \tau = a^\top\Sigma_U b/A$.

*$K$ is symmetric PSD.* From $K = \Sigma_U - \Sigma_U\Sigma_W^{-1}\Sigma_U = \Sigma_U^{1/2}\big(I - \Sigma_U^{1/2}\Sigma_W^{-1}\Sigma_U^{1/2}\big)\Sigma_U^{1/2}$ and $\Sigma_W \succeq \Sigma_U$ (so $\Sigma_U^{1/2}\Sigma_W^{-1}\Sigma_U^{1/2}\preceq I$), the middle factor is PSD, hence $K \succeq 0$. The $\Sigma_U=I$ form (5) is immediate: $K = \Sigma_W^{-1}\Sigma_\delta = \mathrm{diag}\big(\delta_j/(1+\delta_j)\big) = \mathrm{diag}(1-R_j)$.

*Scalar reduction.* At $k=1$ with $\Sigma_U=u$, $\Sigma_\delta=w$: $K = uw/(u+w)$, so $B_W = bauw/(a^2uw+s(u+w))$, matching [8]. $\qquad\blacksquare$

**Theorem 2.** *(b) $\Rightarrow$ (a).* If $\Sigma_\delta = \rho\Sigma_U$ then $\Sigma_W = (1+\rho)\Sigma_U$ and $K = \Sigma_U\Sigma_W^{-1}\Sigma_\delta = \tfrac{\rho}{1+\rho}\Sigma_U =: \lambda\Sigma_U$. So $B_W = \lambda\,a^\top\Sigma_U b/(s+\lambda a^\top\Sigma_U a)$, sharing the sign of $B_{\text{naive}}$, and
$$
\frac{B_W}{B_{\text{naive}}} = \frac{\lambda\,(s + \|a\|^2_{\Sigma_U})}{s + \lambda\|a\|^2_{\Sigma_U}} = \frac{\lambda\|a\|^2_{\Sigma_U} + \lambda s}{\lambda\|a\|^2_{\Sigma_U} + s} \le 1
$$
since $\lambda \le 1$. Thus (a) holds.

*(a) $\Rightarrow$ (b).* Assume (a). Fixing $a,b$ and letting $s\to\infty$, both biases behave like (numerator)$/s$, so $|B_W|\le|B_{\text{naive}}|$ for all $s$ forces
$$
|a^\top K b| \;\le\; |a^\top \Sigma_U b| \qquad\text{for all } a,b\in\mathbb R^k.
\tag{6}
$$
In particular $a^\top\Sigma_U b = 0 \Rightarrow a^\top K b = 0$. Fix $b\ne0$: the linear functional $a\mapsto a^\top(Kb)$ vanishes on the hyperplane $\{a: a^\top(\Sigma_U b)=0\}$, so $Kb$ is parallel to $\Sigma_U b$, i.e. $Kb = \mu(b)\,\Sigma_U b$. Writing $v = \Sigma_U b$ (which ranges over all of $\mathbb R^k$ as $b$ does), $K\Sigma_U^{-1}v = \mu\,v$ for every $v$, so $K\Sigma_U^{-1} = \lambda I$ for a single scalar $\lambda$, i.e. $K = \lambda\Sigma_U$. Feeding $K=\lambda\Sigma_U$ back into (6) gives $|\lambda|\le1$, and $K\succeq0$, $\Sigma_U\succ0$ give $\lambda\ge0$. Finally $K=\lambda\Sigma_U$ means $\Sigma_U\Sigma_W^{-1}\Sigma_\delta = \lambda\Sigma_U \Rightarrow \Sigma_W^{-1}\Sigma_\delta = \lambda I \Rightarrow \Sigma_\delta = \lambda\Sigma_W = \lambda(\Sigma_U+\Sigma_\delta) \Rightarrow \Sigma_\delta = \tfrac{\lambda}{1-\lambda}\Sigma_U$, which is (b) with $\rho=\lambda/(1-\lambda)$. $\qquad\blacksquare$

**Corollary 3.** Work in the $\Sigma_U=I$ case (already a witness), so $B_{\text{naive}} = \big(\textstyle\sum_j a_jb_j\big)/(s+\|a\|^2)$ and, by (5), $B_W = \big(\sum_j a_jb_j(1-R_j)\big)/(s+\sum_j a_j^2(1-R_j))$.

*(iii) Manufacture.* Take $k=2$, $a=(1,1)$, $b=(1,-1)$, so $\sum_j a_jb_j = 0$ and $B_{\text{naive}}=0$. Then $a^\top K b = (1-R_1)-(1-R_2) = R_2 - R_1$, nonzero iff $R_1\ne R_2$. With $R_1=0.9,R_2=0.3,s=1$: $B_W = (0.1-0.7)/(1+0.1+0.7) = -0.6/1.8 = -1/3$. The oracle and naive are both unbiased; adjustment alone produces bias $-1/3$.

*(ii) Reversal.* Perturb to $b=(1,-\tfrac12)$: $B_{\text{naive}} = \tfrac12/(s+2) > 0$, while $a^\top K b = (1-R_1) - \tfrac12(1-R_2)$, which is negative once $R_2$ is small enough (e.g. $R_1=0.9,R_2=0.2$: $0.1 - 0.4 = -0.3<0$). Then $B_W<0<B_{\text{naive}}$.

*(i) Amplification, unbounded ratio.* Let the biases nearly cancel, $a=(2,1.9)$, $b=(1,-1)$, $s=1$, so $B_{\text{naive}} = (2-1.9)/(1+4+3.61)=0.1/8.61\approx0.0116$. Make the strongly biasing confounder poorly measured and the offsetting one well measured, $R_1=0.1,R_2=0.9$: $a^\top K b = 2(0.9) - 1.9(0.1) = 1.8-0.19 = 1.61$, and $B_W = 1.61/(1+4(0.9)+3.61(0.1)) = 1.61/4.961\approx0.3245$, so $B_W/B_{\text{naive}}\approx 28$. Sending the residual naive bias to $0$ (with $a^\top Kb$ held fixed) sends the ratio to $\infty$. $\qquad\blacksquare$

**The mechanism in one line.** Adjustment shrinks confounder $j$'s bias contribution by its reliability, keeping only the fraction $1-R_j$. The naive estimand keeps all of every contribution ($1-R_j\equiv1$) and the oracle keeps none ($1-R_j\equiv0$); both weight the confounders *uniformly*, so a cancellation $\sum_j a_jb_j=0$ that holds for one holds for the other. The proxy weights them *non*-uniformly, by $1-R_j$, and a uniform cancellation is not a non-uniform one. Homogeneous reliability ($R_j\equiv R$) is exactly the restored-uniformity case (b).

## 4. Experiments

`sim.py` (seed `20260702`) verifies (4) two ways and exhibits the failure modes. First, it draws $2\times10^5$ random parameter sets (dimension $k=1,\dots,4$, random $\Sigma_U\succ0$, diagonal $\Sigma_\delta$, random $a,b,s$) and checks the closed form (4) against the direct population block-inverse estimand: the maximum discrepancy is $1.5\times10^{-14}$. Second, it scans $3\times10^5$ random $(a,b,s)$ under each reliability regime and records the worst amplification $|B_W|-|B_{\text{naive}}|$:

| reliability regime | worst $|B_W|-|B_{\text{naive}}|$ | any sign reversal? |
|---|---|---|
| **homogeneous** ($\Sigma_\delta \propto \Sigma_U$) | $-9.8\times10^{-8}$ (i.e. none) | no |
| **heterogeneous** (diagonal $\Sigma_\delta$) | $+2.36$ | yes |

Homogeneous reliability never amplifies and never reverses, to numerical tolerance — Theorem 2 in action. Heterogeneous reliability does both. The named scenarios, with population theory and finite-sample OLS ($n=2\times10^5$, $40$ replications) side by side:

| scenario | $B_{\text{naive}}$ (thy / MC) | $B_W$ (thy / MC) | $B_{\text{oracle}}$ (thy / MC) |
|---|---|---|---|
| scalar $k=1$, $R=0.5$ | $+0.500 / +0.500$ | $+0.333 / +0.333$ | $0 / -0.000$ |
| $k=2$ homogeneous ($\Sigma_\delta=\Sigma_U$) | $+0.217 / +0.217$ | $+0.170 / +0.170$ | $0 / +0.000$ |
| $k=2$ **manufacture** ($R{=}.9,.3$) | $+0.000 / +0.000$ | $\mathbf{-0.333} / -0.334$ | $0 / +0.000$ |
| $k=2$ **reversal** ($R{=}.9,.2$) | $+0.167 / +0.167$ | $\mathbf{-0.158} / -0.158$ | $0 / -0.001$ |
| $k=2$ **amplify** ($R{=}.1,.9$) | $+0.012 / +0.012$ | $\mathbf{+0.325} / +0.325$ | $0 / +0.000$ |

Every theoretical bias matches Monte Carlo to three digits, and every oracle bias is zero to Monte-Carlo error. The scalar and homogeneous rows are safe (proxy strictly between naive and oracle); the three bold rows are the impossibilities of the scalar world: bias created from nothing, bias flipped in sign, and bias magnified $28$-fold.

![Two uncorrelated confounders whose bias contributions exactly cancel, so the **naive estimator is unbiased everywhere** ($B_{\text{naive}}\equiv0$). Color is the proxy-adjusted bias $B_W$ over the reliability grid $(R_1,R_2)$. Adjustment is harmless only on the diagonal $R_1=R_2$ (dashed) — homogeneous reliability, Theorem 2. Off it, adjusting for the two mismeasured confounders manufactures bias, of a sign set by which proxy is the more reliable.](figs/manufacture_heatmap.png)

![Confounder 1 held well measured ($R_1=0.9$); confounder 2's reliability $R_2$ swept downward. The scalar theory confines $B_W/B_{\text{naive}}$ to the green band $(0,1)$ for all reliabilities. With a second confounder the curve leaves the band: it crosses $0$ (sign reversal) and drops below $-1$ (amplification, $|B_W|>|B_{\text{naive}}|$).](figs/overshoot_ratio.png)

## 5. Discussion

The scalar note this sequel builds on ended its limitations with a promissory note: "with vector confounders the scalar $R_X^2$ becomes a matrix object and the bound generalizes to a partial-$R^2$ statement that I have not worked out here." Working it out reverses the headline. For one confounder, reliability is a *guarantee*: whatever you cannot measure well, adjusting for your noisy proxy still moves you toward the truth, capped by $R$. For two or more, reliability is only a guarantee when it is **uniform across the confounders** (relative to their covariance); heterogeneous reliability — the generic case, since different variables are measured by different instruments — voids it. The exact statement (4) shows why: adjustment does not attenuate the *aggregate* bias, it attenuates each confounder's contribution *separately*, by that confounder's own reliability, and separate attenuation is not the same as joint attenuation once the contributions can offset one another.

That adjusting for a mismeasured confounder can increase or reverse bias is not itself new: it is known qualitatively for a single non-differentially mismeasured *ordinal* confounder (Ogburn & VanderWeele [3]), for the interplay of exposure and confounder reliability (Zoh et al. [5], who trace an apparent non-monotonicity to an equal-reliability assumption), and in simulation studies of residual confounding [6] and of adjustment sets in machine learning [7]. The general multivariate errors-in-variables reliability-matrix machinery is Gleser's [1]. What is new here is the **exact closed form for the partialling-out estimand with an arbitrary reliability structure** (4)–(5), and the **sharp characterization of when safety holds** (Theorem 2): not a sign condition on a particular configuration, but an *iff* on the error geometry that quantifies over all effect vectors at once, together with the clean separation of the three failure modes and the observation — absent from the scalar picture — that adjustment can be strictly worse than doing nothing (rows 3–5). The $\Sigma_U=I$ reliability-weighted form (5) makes the diagnostic operational: bias survives adjustment in proportion to each confounder's *un*reliability $1-R_j$, so a single badly-measured confounder can dominate the residual even when everything else is pristine.

Practically, the message is a caution against the reflex "I controlled for it, so it's handled" when *it* is several proxies of differing quality. If reliabilities are known (from validation substudies), (4) is directly computable and turns into an honest residual-bias number and a sign check; if only their spread is known, homogeneity is the property to argue for before trusting the adjustment — and homogeneity, not mere high reliability, is what buys back the scalar guarantee.

**Limitations.** The clean formulas are, as before, creatures of the linear–Gaussian world with additive effects and no $U\times X$ interaction. (i) I take a **scalar treatment** $X$ measured without error; error in $X$ adds the classical attenuation of the *exposure* on top of this and is not modeled. (ii) Error is **classical**: additive, mean-zero, independent, and coordinate-wise diagonal $\Sigma_\delta$; correlated measurement errors (non-diagonal $\Sigma_\delta$) are still covered by the matrix form (4) but not by the reliability-weighted reading (5), and differential error (error correlated with $Y$) breaks the setup entirely. (iii) The results are about **population bias**, not variance: adjusting for near-noise proxies also inflates the sampling variance of $\hat\tau$, a second cost the bias accounting is silent on, so "manufacture" of a small bias may still be worthwhile if it removes a large one elsewhere — the note quantifies the bias ledger, not the mean-squared-error trade. (iv) $\Sigma_U\succ0$ (no perfectly collinear confounders) is used throughout; the rank-deficient case is Gleser's non-identifiability [1]. (v) The dichotomy is about *guaranteed* safety over all $(a,b,s)$; for a specific known configuration, heterogeneous reliability may still happen to be safe — the theorem says only that nothing protects you in general.

## References

1. L. J. Gleser. "The importance of assessing measurement reliability in multivariate regression." *Journal of the American Statistical Association*, 87(419):696–707, 1992. doi:10.1080/01621459.1992.10475271.
2. R. J. Carroll, D. Ruppert, L. A. Stefanski, and C. M. Crainiceanu. *Measurement Error in Nonlinear Models: A Modern Perspective*, 2nd ed. Chapman & Hall/CRC, 2006.
3. E. L. Ogburn and T. J. VanderWeele. "Bias attenuation results for nondifferentially mismeasured ordinal and coarsened confounders." *Biometrika*, 100(1):241–248, 2013. doi:10.1093/biomet/ass054.
4. A general condition for bias attenuation by a nondifferentially mismeasured confounder. *Biometrika*, 112(3):asaf026, 2025. doi:10.1093/biomet/asaf026 (preprint arXiv:2409.12928).
5. R. S. Zoh, D. M. Thomas, C. D. Tekwe, X. Yu, C. J. Vorland, N. V. Dhurandhar, D. M. Klurfeld, and D. B. Allison. "Adjusting for covariates representing potential confounders, mediators, or competing predictors in the presence of measurement error: dispelling a potential misapprehension and insights for optimal study design." *F1000Research*, 13:827, 2024/2025. doi:10.12688/f1000research.152466.2.
6. C. L. Fewell, G. Davey Smith, and J. A. C. Sterne. "The impact of residual and unmeasured confounding in epidemiologic studies: a simulation study." *American Journal of Epidemiology*, 166(6):646–655, 2007. doi:10.1093/aje/kwm165.
7. K. Sinha and P. R. Raamana. "Indiscriminate adjustment for confounders is worse than you think and what can be done about it." *Alzheimer's & Dementia*, 20(S2):e090913, 2024. doi:10.1002/alz.090913.
8. E. Nakamura-Lindqvist. "Adjusting for a noisy confounder removes at most a fraction $R$ of the bias — and usually less." *Chase's Journal*, June 2026 (`submissions/2026-06/m7-proxy-confounder-reliability-bound`).
9. R. Frisch and F. V. Waugh. "Partial time regressions as compared with individual trends." *Econometrica*, 1(4):387–401, 1933. (The Frisch–Waugh–Lovell partialling-out theorem; M. C. Lovell, *JASA*, 58:993–1010, 1963.)

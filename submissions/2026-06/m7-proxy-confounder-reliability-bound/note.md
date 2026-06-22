# Adjusting for a noisy confounder removes at most a fraction $R$ of the bias — and usually less

**Author.** Esther Nakamura-Lindqvist — a (fictional) researcher in causal inference and measurement error who distrusts any "adjusted" estimate whose control variables were measured with a ruler.

**Submitted to *Chase's Journal*.** 2026-06-22

## Abstract

When a confounder cannot be measured exactly, the reflex is to adjust for a noisy proxy and hope that most of the bias goes away. *How much* goes away? In the linear–Gaussian model with a classical-error proxy of reliability $R = \mathrm{Var}(U)/\mathrm{Var}(W)$, I give an exact answer: the fraction of confounding bias removed by adjusting for the proxy is
$$ F \;=\; \frac{R\,(1-R_X^2)}{R\,(1-R_X^2) + (1-R)}, \tag{1}$$
where $R_X^2$ is the share of the treatment's variance explained by the confounder. The clean consequence is a sharp bound: $F \le R$, with equality **iff** the confounder has no effect on treatment ($R_X^2 = 0$). So a proxy that is, say, $90\%$ reliable removes *strictly less* than $90\%$ of the bias whenever it is genuinely confounding, and the shortfall grows with how strongly the confounder drives treatment. The result sharpens the qualitative "adjustment attenuates bias" literature into a quantitative rule of thumb and a worst-case reading, and a simulation confirms every formula to Monte-Carlo precision.

## 1. Setup and notation

A single unmeasured confounder $U$ sits behind a treatment $X$ and an outcome $Y$. We want the causal effect $\tau$ of $X$ on $Y$. The data-generating process is the standard linear–Gaussian structural model with mean-zero, mutually independent noise terms:
$$
\begin{aligned}
U &\sim \mathcal N(0,\,u),\\
X &= a\,U + \varepsilon_X, &\varepsilon_X &\sim \mathcal N(0,\,s),\\
Y &= \tau\,X + b\,U + \varepsilon_Y, &\varepsilon_Y &\sim \mathcal N(0,\,v).
\end{aligned}
\tag{2}
$$
Here $\tau$ is the target (the structural coefficient of $X$), and $b \ne 0$, $a \ne 0$ make $U$ a genuine confounder: it affects both $X$ and $Y$. We cannot observe $U$. Instead we observe a **classical-error proxy**
$$ W = U + \varepsilon_W, \qquad \varepsilon_W \sim \mathcal N(0,\,w), \tag{3}$$
with $\varepsilon_W$ independent of everything else (nondifferential, uncorrelated error). Two summaries of the model will carry the whole story:
$$ R \;=\; \frac{\mathrm{Var}(U)}{\mathrm{Var}(W)} \;=\; \frac{u}{u+w} \in (0,1], \qquad R_X^2 \;=\; \frac{\mathrm{Var}(aU)}{\mathrm{Var}(X)} \;=\; \frac{a^2 u}{a^2 u + s} \in [0,1). \tag{4}$$
$R$ is the **reliability** of the proxy — the fraction of its variance that is signal. $R_X^2$ is the population $R^2$ of regressing $X$ on $U$: the fraction of the treatment's variance the confounder explains. Write $A := \mathrm{Var}(X) = a^2 u + s$.

We compare three population ordinary-least-squares estimands of the coefficient on $X$:

- **naive**: regress $Y$ on $X$ alone;
- **proxy**: regress $Y$ on $(X, W)$ — adjust for what we can measure;
- **oracle**: regress $Y$ on $(X, U)$ — adjust for the truth.

The oracle recovers $\tau$ exactly (it blocks the back-door path $X \leftarrow U \rightarrow Y$). The question is where the proxy lands between the naive estimand and the oracle.

## 2. Contribution

**Proposition (exact residual confounding).** *Under model (2)–(3), the population OLS coefficient on $X$ has bias*
$$
\underbrace{B_{\mathrm{full}} \;=\; \frac{b\,a\,u}{A}}_{\text{naive, no adjustment}}
\qquad\text{and}\qquad
\underbrace{B_W \;=\; \frac{b\,a\,u\,w}{a^2 u\,w + s\,(u+w)}}_{\text{adjusting for the proxy }W}. \tag{5}
$$
*Consequently the fraction of the naive bias removed by adjusting for $W$ is*
$$
F \;:=\; 1 - \frac{B_W}{B_{\mathrm{full}}} \;=\; \frac{s\,u}{s\,u + w\,A} \;=\; \frac{R\,(1-R_X^2)}{R\,(1-R_X^2) + (1-R)}. \tag{6}
$$

**Corollary (the reliability bound).** *$\,0 < F \le R$, with equality iff $R_X^2 = 0$. Moreover $F$ is strictly increasing in $R$ and strictly decreasing in $R_X^2$.*

The corollary is the take-home message. A practitioner who knows their proxy is "$90\%$ reliable" instinctively expects to remove $90\%$ of the bias. Equation (6) says they remove $90\%$ only in the degenerate case where $U$ does not affect treatment at all (and then there was nothing to remove via the $X$-path anyway). The instant the confounder actually moves treatment ($R_X^2 > 0$), the proxy underperforms its reliability, and the shortfall $R - F$ widens as $R_X^2 \to 1$.

## 3. Argument

**Bias of the naive estimand.** With one regressor, the OLS coefficient on $X$ is $\mathrm{Cov}(Y,X)/\mathrm{Var}(X)$. From (2),
$$ \mathrm{Cov}(Y,X) = \tau A + b\,\mathrm{Cov}(U,X) = \tau A + b\,a u, $$
so the coefficient is $\tau + b a u / A$, giving $B_{\mathrm{full}} = b a u / A$ — the textbook omitted-variable bias.

**Bias of the proxy-adjusted estimand.** Use the Frisch–Waugh–Lovell theorem: the coefficient on $X$ in the regression of $Y$ on $(X,W)$ equals $\mathrm{Cov}(Y,\tilde X)/\mathrm{Var}(\tilde X)$, where $\tilde X$ is the residual of $X$ after projecting on $W$. The projection coefficient is $\mathrm{Cov}(X,W)/\mathrm{Var}(W) = a u/(u+w)$ (since $\mathrm{Cov}(X,W) = \mathrm{Cov}(aU+\varepsilon_X,\,U+\varepsilon_W) = a u$). Hence
$$
\mathrm{Var}(\tilde X) = A - \frac{(a u)^2}{u+w} = a^2 u\Big(1 - \tfrac{u}{u+w}\Big) + s = \frac{a^2 u\, w}{u+w} + s. \tag{7}
$$
For the numerator, $\mathrm{Cov}(Y,W) = \tau\,\mathrm{Cov}(X,W) + b\,\mathrm{Cov}(U,W) = \tau a u + b u$, so
$$
\mathrm{Cov}(Y,\tilde X) = \mathrm{Cov}(Y,X) - \frac{\mathrm{Cov}(X,W)}{\mathrm{Var}(W)}\mathrm{Cov}(Y,W)
= \big(\tau A + bau\big) - \frac{au}{u+w}\big(\tau a u + b u\big).
$$
Collecting the $\tau$-terms reproduces $\tau\,\mathrm{Var}(\tilde X)$ exactly (as it must), and the $b$-terms give
$$
b a u - \frac{a u}{u+w}\,b u = b a u\Big(1 - \tfrac{u}{u+w}\Big) = \frac{b a u\, w}{u+w}. \tag{8}
$$
Dividing (8) by (7),
$$
B_W = \frac{b a u\, w/(u+w)}{a^2 u w/(u+w) + s} = \frac{b a u\, w}{a^2 u w + s(u+w)},
$$
which is (5). Two sanity checks: $w \to 0$ (perfect proxy, $W=U$) gives $B_W \to 0$; $w \to \infty$ (pure noise) gives $B_W \to bau/(a^2u+s) = B_{\mathrm{full}}$. Both are correct.

**The fraction removed.** Dividing,
$$
\frac{B_W}{B_{\mathrm{full}}} = \frac{b a u w}{a^2 uw + s(u+w)}\cdot\frac{A}{bau}
= \frac{w A}{a^2 u w + s(u+w)} = \frac{wA}{w(a^2u+s) + su} = \frac{wA}{wA + su},
$$
so $F = 1 - B_W/B_{\mathrm{full}} = su/(su + wA)$. Substituting $w = u(1-R)/R$, $s = A(1-R_X^2)$, and $a^2u = A R_X^2$ turns this into the reliability form (6); the algebra is a one-liner and is checked symbolically against the raw form over $10^5$ random parameter vectors in `sim.py` (agreement to $4\times10^{-15}$).

**Proof of the bound.** Write $F = 1/\big(1 + \tfrac{wA}{su}\big)$ and note $A/s = (a^2u+s)/s = 1/(1-R_X^2) \ge 1$. Then
$$
\frac{wA}{su} = \frac{w}{u}\cdot\frac{A}{s} = \frac{1-R}{R}\cdot\frac{1}{1-R_X^2} \;\ge\; \frac{1-R}{R},
$$
with equality iff $R_X^2 = 0$. Since $F$ is decreasing in this quantity,
$$
F \;\le\; \frac{1}{1 + (1-R)/R} = R,
$$
with equality iff $R_X^2=0$. Monotonicity in $R$ and in $R_X^2$ is immediate from the same display. $\qquad\blacksquare$

**Why $R_X^2$ and not just reliability.** There is a clean interpretation. Adjusting for $W$ is, in the population, equivalent to adjusting for $\hat U := \mathbb E[U \mid W] = R\,W$, the best linear predictor of the confounder. What remains unblocked is the residual confounder $U^\perp = U - \hat U$, with variance $\mathrm{Var}(U^\perp) = u(1-R)$. So a reliability-$R$ proxy leaves behind an unmeasured confounder carrying a $(1-R)$ share of the original confounder's variance. The damage this residual does is **amplified** by how tightly $U$ is coupled to $X$: when $U$ explains most of the treatment's variation ($R_X^2$ near $1$), $\tilde X$ — the part of treatment orthogonal to the measured proxy — is almost pure leftover confounder, so the residual confounding is barely dented. That coupling is exactly the factor $1/(1-R_X^2) = A/s$ in the proof. Reliability tells you how much confounder *variance* you captured; $R_X^2$ tells you how much that capture is *worth* for de-confounding treatment.

## 4. Experiments

`sim.py` (seed `20260622`) does two things. First, it checks (5)–(6) against Monte-Carlo OLS: for each of three scenarios it draws $n=40{,}000$ observations, $200$ times, and runs the naive, proxy, and oracle regressions. Second, it draws $10^5$ random parameter vectors to verify the identity (6) and the bound $F\le R$ globally.

| scenario | $R$ | $R_X^2$ | $B_{\mathrm{full}}$ (thy / MC) | $B_W$ (thy / MC) | $F$ (thy / MC) |
|---|---|---|---|---|---|
| moderate | 0.667 | 0.500 | 0.500 / 0.500 | 0.250 / 0.250 | 0.500 / 0.500 |
| $U$ drives $X$ | 0.500 | 0.941 | 0.706 / 0.706 | 0.667 / 0.667 | **0.056** / 0.055 |
| $U$ weakly drives $X$ | 0.870 | 0.143 | 0.286 / 0.285 | 0.043 / 0.043 | 0.851 / 0.851 |

The oracle bias is $0$ to three Monte-Carlo standard errors in every scenario, as it must be. The middle row is the cautionary tale: a proxy of reliability $R = 0.5$ — half signal, half noise, not obviously useless — removes only **$5.6\%$** of the confounding bias, because the confounder explains $94\%$ of the treatment's variance. Over the $10^5$-vector grid, the identity (6) holds to $4\times10^{-15}$ and $F \le R$ holds in $100\%$ of draws (largest violation $-2.7\times10^{-10}$, i.e. none).

![Fraction of confounding bias removed, $F$, versus proxy reliability $R$, for several values of $R_X^2$. The dashed line $F=R$ is the naive expectation; every curve lies on or below it, touching only at $R_X^2=0$. Open circles are the Monte-Carlo scenarios from the table.](figs/fraction_removed.png)

![The shortfall $R-F$ — bias you expected your proxy to remove but didn't — as a function of reliability $R$ and the treatment's confounder-share $R_X^2$. White contours mark fixed values of $F$. The bottom edge ($R_X^2=0$) is the only place the proxy delivers its full reliability; performance degrades sharply as the confounder takes over the treatment.](figs/shortfall_heatmap.png)

## 5. Discussion

The qualitative fact that adjusting for a nondifferentially mismeasured confounder *attenuates* bias toward the truth (rather than overshooting or reversing it) has a careful literature, most of it phrased as sign or monotonicity conditions — Ogburn & VanderWeele [3] for ordinal/coarsened confounders, and more recently a general attenuation condition for additive-noise and coarsened proxies [4]. Equation (6) is consistent with that picture — $F \in (0,1)$ means the proxy estimand always lands strictly between naive and oracle — but it answers a sharper question: not *whether* bias attenuates, but *by what fraction*, in closed form. The bound $F \le R$ packages this into a one-number worst case that needs only the proxy's reliability, a quantity that reliability/validation studies routinely report.

The mechanism is the same one behind classical attenuation (regression dilution) of a mismeasured *exposure* [1,2], but the accounting is different and, I think, less widely appreciated: for a mismeasured *control*, the relevant amplifier is the confounder's grip on treatment, $R_X^2$, not just the proxy's reliability. This reframes a common modeling instinct. Faced with a confounder measured with error, analysts often reason "the proxy is decent, so most of the confounding is handled." Equation (6) says the right diagnostic pairs reliability with $R_X^2$: a highly reliable proxy of a confounder that only weakly perturbs treatment is nearly sufficient (top of Fig. 1), whereas even a good proxy of a confounder that dominates treatment selection is close to useless (middle row of the table). This also gives a quick sensitivity check — if you can bound the reliability $R$ from a validation substudy, then $R$ is a hard ceiling on the bias fraction you have removed, and the leftover $(1-F)\,B_{\mathrm{full}}$ can be reported as residual confounding.

**Limitations.** The clean formula is a creature of the linear–Gaussian world. It assumes (i) linear structural equations with additive effects and no $U\times X$ interaction; (ii) a single scalar confounder and a single scalar proxy; (iii) **classical** measurement error — additive, mean-zero, and independent of $U$, $X$, $Y$, and all noise (nondifferential and uncorrelated). Departures matter: differential error (proxy error correlated with the outcome) can *reverse* the sign of the adjustment and break $F \le R$ entirely; Berkson error obeys different algebra; with vector confounders the scalar $R_X^2$ becomes a matrix object and the bound generalizes to a partial-$R^2$ statement that I have not worked out here. The result is about *bias* of the population estimand, not finite-sample variance — adjusting for a near-noise proxy also inflates variance, a second cost not modeled here. Finally, $F \le R$ is reassuring as a ceiling but unhelpful as a floor: the actual $F$ can be far below $R$, so reliability alone never licenses the claim "most of the confounding is handled" without also pinning down $R_X^2$.

## References

1. J. M. Frost and S. G. Thompson. "Correcting for regression dilution bias: comparison of methods for a single predictor variable." *Journal of the Royal Statistical Society: Series A*, 163(2):173–189, 2000. doi:10.1111/1467-985X.00164.
2. R. J. Carroll, D. Ruppert, L. A. Stefanski, and C. M. Crainiceanu. *Measurement Error in Nonlinear Models: A Modern Perspective*, 2nd ed. Chapman & Hall/CRC, 2006.
3. E. L. Ogburn and T. J. VanderWeele. "Bias attenuation results for nondifferentially mismeasured ordinal and coarsened confounders." *Biometrika*, 100(1):241–248, 2013. doi:10.1093/biomet/ass054.
4. A general condition for bias attenuation by a nondifferentially mismeasured confounder. *Biometrika*, 112(3):asaf026, 2025. doi:10.1093/biomet/asaf026 (preprint arXiv:2409.12928).
5. C. L. Fewell, G. Davey Smith, and J. A. C. Sterne. "The impact of residual and unmeasured confounding in epidemiologic studies: a simulation study." *American Journal of Epidemiology*, 166(6):646–655, 2007. doi:10.1093/aje/kwm165.
6. J. Pearl. *Causality: Models, Reasoning, and Inference*, 2nd ed. Cambridge University Press, 2009.

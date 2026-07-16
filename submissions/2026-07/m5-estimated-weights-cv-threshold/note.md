# When not to weight: a sharp estimation penalty for plug-in inverse-variance weighting, and a coefficient-of-variation threshold for just averaging

**Author.** Elias Kandelaars — applied probabilist working on estimation and meta-analysis; keeps a soft spot for the estimators nobody bothers to prove things about.
**Submitted to *Chase's Journal*.** 2026-07-16

## Abstract

Inverse-variance weighting (IVW) is the reflex for combining $k$ independent unbiased estimates of a common quantity: weight each by its precision and you get the minimum-variance linear combination. But the precisions are almost never known — they are estimated from the same data — and the resulting *plug-in* estimator (the Graybill–Deal estimator) is strictly worse than it looks. We give a clean, exact statement of the price. When the true variances are equal, plug-in IVW is unbiased but its variance exceeds that of a plain average by the factor $1 + \tfrac{2(k-1)}{k(n-1)} + O(n^{-2})$, where $n$ is the per-group sample size; the leading term is a genuine asymptotic coefficient, not a bound. Turning this into a decision rule, plug-in IVW beats simple averaging only when the true variances are heterogeneous enough — precisely, when their squared coefficient of variation exceeds the same threshold, $c_v^2 \gtrsim \tfrac{2(k-1)}{k(n-1)}$. A simulation confirms the penalty constant to three digits across $k$ and $n$ and locates the crossover within a few percent of the predicted threshold. The upshot is a memorable rule of thumb: **with $n$ observations per group, don't estimate weights unless the variances differ by more than a squared CV of about $2/n$ — otherwise just average.**

## 1. Setup and notation

We observe $k$ independent groups. Group $i$ contributes $n$ i.i.d. observations
$X_{i1},\dots,X_{in} \sim \mathcal N(\theta, \sigma_i^2)$ with a **common mean**
$\theta$ (the estimand) and its own unknown variance $\sigma_i^2$. Write the group
mean and unbiased sample variance
$$
\bar X_i = \frac1n\sum_{j=1}^n X_{ij}, \qquad s_i^2 = \frac1{n-1}\sum_{j=1}^n (X_{ij}-\bar X_i)^2 .
$$
Each $\bar X_i$ is an unbiased estimate of $\theta$ with variance $\sigma_i^2/n$. We
compare three ways of combining them into a single estimate $\hat\theta = \sum_i w_i \bar X_i$ with $\sum_i w_i = 1$:

- **EQUAL** — the plain average, $w_i = 1/k$.
- **ORACLE** — inverse-variance weights with *known* variances, $w_i^\star \propto 1/\sigma_i^2$. This is the minimum-variance unbiased linear combination (Gauss–Markov / Aitken), with variance $\big(\sum_i n/\sigma_i^2\big)^{-1}$.
- **PLUGIN** — feasible IVW with *estimated* variances, $\hat w_i \propto 1/s_i^2$. For $k$ populations this is exactly the **Graybill–Deal estimator** [1] of the common mean.

Everything below is about the gap between PLUGIN and its two idealized neighbors. Let $\nu = n-1$ be the per-group degrees of freedom, and let
$$
\overline{\sigma^2} = \frac1k\sum_i \sigma_i^2, \qquad
c_v^2 = \frac{\operatorname{Var}_i(\sigma_i^2)}{\overline{\sigma^2}^{\,2}} = \frac{\frac1k\sum_i(\sigma_i^2 - \overline{\sigma^2})^2}{\overline{\sigma^2}^{\,2}}
$$
denote the (squared) **coefficient of variation of the true variances** — a scale-free measure of how heteroscedastic the groups are. Homoscedastic means $c_v^2 = 0$.

## 2. The homoscedastic penalty

Two structural facts make the Gaussian case exactly analyzable. First, for Gaussian data $\bar X_i \perp s_i^2$ (Basu / Cochran), so the plug-in *weights are independent of the group means they multiply*. Second, the weights depend on the variances only through ratios, so under homoscedasticity ($\sigma_i^2 \equiv \sigma^2$) they carry no information about which group is more precise — they are pure noise.

**Proposition 1 (unbiased, but penalized).**
*The plug-in estimator is exactly unbiased, $\mathbb E[\hat\theta_{\rm PLUGIN}] = \theta$, and*
$$
\operatorname{Var}(\hat\theta_{\rm PLUGIN}) = \mathbb E\!\Big[\textstyle\sum_i \hat w_i^2\, \sigma_i^2/n\Big]. \tag{1}
$$
*In the homoscedastic case this is $(\sigma^2/n)\,\mathbb E[\sum_i \hat w_i^2]$, and since $\sum_i \hat w_i^2 \ge 1/k$ with equality iff all weights are equal,*
$$
\frac{\operatorname{Var}(\hat\theta_{\rm PLUGIN})}{\operatorname{Var}(\hat\theta_{\rm EQUAL})} = k\,\mathbb E\!\Big[\textstyle\sum_i \hat w_i^2\Big] \;>\; 1, \tag{2}
$$
*strictly, whenever the weights are genuinely random (any finite $n$).*

*Proof.* Unbiasedness: $\hat w = (\hat w_i)$ is a function of $(s_1^2,\dots,s_k^2)$ alone, hence independent of $(\bar X_1,\dots,\bar X_k)$, and $\sum_i \hat w_i = 1$, so $\mathbb E[\hat\theta \mid \hat w] = \sum_i \hat w_i\, \theta = \theta$. For the variance, condition on $\hat w$: the $\bar X_i$ are independent with variance $\sigma_i^2/n$, so $\operatorname{Var}(\hat\theta\mid\hat w) = \sum_i \hat w_i^2\,\sigma_i^2/n$. Since $\mathbb E[\hat\theta\mid\hat w]=\theta$ is constant in $\hat w$, the law of total variance gives (1). The homoscedastic specialization uses $\sigma_i^2\equiv\sigma^2$; the inequality $\sum_i \hat w_i^2 \ge (\sum_i \hat w_i)^2/k = 1/k$ is Cauchy–Schwarz (equivalently, power-mean), with equality iff the $\hat w_i$ are all equal, an event of probability zero for continuous $s_i^2$. The EQUAL estimator has variance $\sigma^2/(kn)$, giving the ratio (2). $\qquad\blacksquare$

So under homoscedasticity — the exact regime where EQUAL is the oracle-optimal choice — estimating weights can only hurt. The next result says *how much*.

**Theorem 2 (the $2(k-1)/(k(n-1))$ penalty).**
*Under homoscedasticity, as $n\to\infty$ with $k$ fixed,*
$$
k\,\mathbb E\!\Big[\textstyle\sum_i \hat w_i^2\Big] \;=\; 1 + \frac{2(k-1)}{k(n-1)} \;+\; O\!\big(n^{-2}\big). \tag{3}
$$
*Equivalently, the relative excess variance of plug-in IVW over simple averaging is $\dfrac{2(k-1)}{k(n-1)}$ to leading order — growing in $k$ (from $\tfrac{1}{n-1}$ at $k=2$ toward $\tfrac{2}{n-1}$ as $k\to\infty$) and shrinking like $1/n$.*

*Proof (leading-order delta method).* Write $s_i^2 = \sigma^2 V_i/\nu$ with $V_i \sim \chi^2_\nu$ i.i.d., so the weights $\hat w_i = (1/V_i)/\sum_j (1/V_j)$ do not depend on $\sigma^2$. Put $V_i = \nu(1+Z_i)$, so $Z_i = (V_i-\nu)/\nu$ has $\mathbb E Z_i = 0$ and $\tau^2 := \operatorname{Var}(Z_i) = 2/\nu$. Expanding $b_i := (1+Z_i)^{-1} = 1 - Z_i + Z_i^2 - \cdots$,
$$
\sum_i \hat w_i^2 = \frac{\sum_i b_i^2}{\big(\sum_i b_i\big)^2}, \qquad
b_i^2 = 1 - 2Z_i + 3Z_i^2 + \cdots .
$$
With $s_1 = \sum_i Z_i$ and $s_2 = \sum_i Z_i^2$, the numerator is $k - 2s_1 + 3s_2 + \cdots$ and $\big(\sum_i b_i\big)^2 = (k - s_1 + s_2 + \cdots)^2 = k^2\big(1 - 2s_1/k + 2s_2/k + s_1^2/k^2 + \cdots\big)$, so
$$
\big(\textstyle\sum_i b_i\big)^{-2} = k^{-2}\big(1 + 2s_1/k - 2s_2/k + 3s_1^2/k^2 + \cdots\big).
$$
Multiplying and keeping terms through second order in the $Z_i$, the first-order pieces cancel and
$$
\sum_i \hat w_i^2 = \frac1k\Big(1 + \tfrac1k s_2 - \tfrac1{k^2}s_1^2 + \cdots\Big).
$$
Take expectations: $\mathbb E[s_2] = k\tau^2$ and $\mathbb E[s_1^2] = \operatorname{Var}(s_1) = k\tau^2$ (independence, $\mathbb E Z_i = 0$), so
$$
\mathbb E\!\Big[\textstyle\sum_i \hat w_i^2\Big] = \frac1k\Big(1 + \tau^2 - \frac{\tau^2}{k}\Big) + O(\tau^4)
= \frac1k\Big(1 + \tau^2\,\frac{k-1}{k}\Big) + O(n^{-2}).
$$
Multiplying by $k$ and substituting $\tau^2 = 2/(n-1)$ gives (3). $\qquad\blacksquare$

A remark on rigor. The expansion is a standard delta-method / Taylor argument; the $O(n^{-2})$ remainder is stated as the correct order but is *not* accompanied by a uniform bound on the neglected moments (the tail of $1/V_i$ has all moments only for $\nu>4$, which we assume). We therefore label (3) a **leading-order asymptotic**, and verify the coefficient — including the claimed limit $(n-1)\cdot(\text{ratio}-1)\to 2(k-1)/k$ — directly by simulation in §4. For $k=2$ the penalty is the especially clean $\tfrac{1}{n-1}$.

## 3. The heteroscedastic crossover: a CV threshold

Now let the variances differ. Two competing effects decide EQUAL vs. PLUGIN.

**The prize for weighting.** Comparing the two *oracle-referenced* variances,
$$
\frac{\operatorname{Var}(\hat\theta_{\rm EQUAL})}{\operatorname{Var}(\hat\theta_{\rm ORACLE})} = \Big(\tfrac1k\textstyle\sum_i \sigma_i^2\Big)\Big(\tfrac1k\sum_i \sigma_i^{-2}\Big) = \overline{\sigma^2}\cdot\overline{\sigma^{-2}} =: H \ge 1, \tag{4}
$$
by Cauchy–Schwarz, with $H=1$ iff homoscedastic. For mild heterogeneity, a second-order expansion of the arithmetic-mean $\times$ arithmetic-mean-of-reciprocals gives the clean approximation
$$
H \;=\; 1 + c_v^2 + O(c_v^4). \tag{5}
$$
So knowing the variances buys you a factor $\approx 1+c_v^2$; the more heterogeneous the groups, the bigger the prize.

**The price for estimating them.** But you don't know the variances. To leading order the estimation penalty of §2 persists at small heterogeneity, inflating PLUGIN relative to ORACLE by the factor $1 + \tfrac{2(k-1)}{k(n-1)}$. Combining, PLUGIN beats EQUAL iff the prize exceeds the price:
$$
\underbrace{H}_{\approx\,1+c_v^2} \;>\; 1 + \frac{2(k-1)}{k(n-1)}
\qquad\Longleftrightarrow\qquad
\boxed{\;c_v^2 \;\gtrsim\; \frac{2(k-1)}{k(n-1)}\;}. \tag{6}
$$

**Decision rule.** *Estimate inverse-variance weights only when the groups' true variances have squared coefficient of variation above $\tfrac{2(k-1)}{k(n-1)}$ (about $2/n$ for many groups). Below that, simple averaging has smaller variance — the noise in the estimated weights costs more than the heterogeneity they could exploit.* For $n$ around 10–20 this threshold is substantial: variances must differ by tens of percent (in CV) before weighting pays.

This is the quantitative content of a folklore warning ("don't trust estimated weights in small samples"): the threshold is not vanishingly small, it scales like $1/n$, and it has a closed form.

## 4. Simulation

`sim.py` (seed `20260716`, ~80 s) checks both claims by direct Monte Carlo; every number below is printed by the script.

**Part A — the penalty constant.** For homoscedastic groups we estimate the exact ratio (2) and compare to (3). Agreement is excellent once $n$ is moderate, and the *coefficient* $2(k-1)/k$ is confirmed by the convergence of $(n-1)(\text{ratio}-1)$:

| $k$ | $n=6$ | $n=11$ | $n=21$ | $n=51$ | $n=101$ | predicted at $n{=}21$ |
|---|---|---|---|---|---|---|
| 2  | 1.167 | 1.091 | 1.048 | 1.020 | 1.010 | $1.050$ |
| 5  | 1.424 | 1.190 | 1.088 | 1.033 | 1.016 | $1.080$ |
| 10 | 1.637 | 1.245 | 1.105 | 1.038 | 1.019 | $1.090$ |

Reading the table: at $n=21$ observations per group, plug-in IVW carries a **4.8% / 8.8% / 10.5%** variance penalty over plain averaging for $k=2/5/10$ — for doing strictly nothing, since the groups are identical. The leading-order formula (3) tracks the simulation to within Monte-Carlo error for $n\gtrsim 20$ and is an underestimate at very small $n$ (where higher-order terms bite, more so for large $k$). The convergence check gives $(n-1)(\text{ratio}-1) \to 1.00,\,1.61,\,1.82$ for $k=2,5,10$ as $n$ grows (targets $1.00,\,1.60,\,1.80$), pinning the coefficient.

![Homoscedastic penalty: markers are simulated relative excess variance of plug-in IVW over simple averaging; solid lines are the closed-form leading term $2(k-1)/(k(n-1))$. Log–log axes; the lines have the correct $-1$ slope in $n$ and the simulation approaches them from above.](figs/penalty.png)

**Part B — the crossover.** Fix $k=8$, $n=11$, so the predicted threshold is $2\cdot7/(8\cdot10)=0.175$. Sweeping heteroscedastic variance profiles of increasing $c_v^2$ (mean fixed to 1), we measure the variance of all three estimators:

| $c_v^2$ | $H$ | EQUAL | PLUGIN | ORACLE | PLUGIN $<$ EQUAL? |
|---|---|---|---|---|---|
| 0.10 | 1.107 | 0.01133 | 0.01254 | 0.01025 | no |
| 0.15 | 1.165 | 0.01141 | 0.01199 | 0.00980 | no |
| 0.175 | 1.196 | 0.01138 | 0.01163 | 0.00951 | no |
| 0.20 | 1.227 | 0.01135 | 0.01130 | 0.00926 | **yes** |
| 0.25 | 1.294 | 0.01138 | 0.01071 | 0.00880 | **yes** |

EQUAL's variance is essentially flat in $c_v^2$ (it equals $\overline{\sigma^2}/(kn)$, and $\overline{\sigma^2}$ is held fixed), ORACLE improves steadily, and PLUGIN crosses EQUAL from above at an empirical $c_v^2 \approx 0.196$ — a few percent above the predicted $0.175$. The small gap is expected and honest: the threshold (6) uses two leading-order approximations ($H\approx 1+c_v^2$, which *underestimates* $H$, and the small-heterogeneity persistence of the §2 penalty, which slightly *underestimates* the finite-$n$ penalty at $n=11$), both of which push the true crossover a little higher. The qualitative and near-quantitative picture holds cleanly.

![When to weight ($k=8$, $n=11$): variance of the simple average (flat), plug-in IVW (falling), and oracle IVW (lower bound) as the true variances become more heterogeneous. Plug-in IVW overtakes simple averaging only to the right of the predicted threshold $c_v^2=2(k{-}1)/(k(n{-}1))=0.175$ (dotted line).](figs/crossover.png)

## 5. Discussion

The message is not that inverse-variance weighting is wrong — with *known* precisions it is optimal, and with strongly heterogeneous groups the feasible version still wins comfortably. The message is that the "obviously better" reflex has a quantifiable cost that practitioners routinely ignore: in the homoscedastic regime plug-in IVW is *strictly dominated* by the estimator that ignores the variances entirely, and the domination extends to a nontrivial band of genuine heterogeneity whose width, $c_v^2 \lesssim 2(k-1)/(k(n-1))$, is set by the per-group sample size. Meta-analyses of few-replicate studies, sensor fusion with short calibration windows, and federated/averaging schemes with small per-shard samples all live in exactly the small-$n$ regime where this band is wide.

**Relation to prior work.** The Graybill–Deal estimator [1] and the efficiency of "weighting by the estimated inverse variance" have a long history: Cochran and Carroll's sampling study [2], exact variance and distribution results for the two-population estimator (Nair [3]), and the broad feasible-GLS/WLS literature showing that estimated weights can be *less* efficient than no weighting under homoscedasticity in small samples [4,5]. The parallel folklore in meta-analysis is that small-sample precision weighting is unreliable and can be beaten by fixed (e.g. sample-size) weights [6,7]. What we add is not the phenomenon but its **sharp, decision-ready form**: (i) the exact homoscedastic ratio (2) and its clean leading constant $2(k-1)/(k(n-1))$ for general $k$, and (ii) the reframing of "should I weight?" as a single scale-free comparison — the coefficient of variation of the variances against that same threshold (6). Exact variance formulas for the Graybill–Deal estimator exist [3] but are intricate; the value here is the memorable closed form and the actionable rule.

**Limitations.** (a) The penalty (3) is a *leading-order* asymptotic in $n$: the coefficient is verified but the $O(n^{-2})$ remainder is not rigorously bounded, and at very small $n$ (say $n\le 6$) it materially understates the true penalty, especially for large $k$. (b) Gaussianity is load-bearing twice — for the exact unbiasedness of PLUGIN (via $\bar X_i \perp s_i^2$) and for the $\chi^2$ law of $s_i^2$; with heavy tails or non-normal data both the independence and the constant change. (c) The crossover (6) inherits two small-heterogeneity approximations, so it is a *rule of thumb*, accurate to a few percent in $c_v^2$ at the sample sizes shown, not an exact boundary; the true crossover sits slightly above it. (d) We assume a genuinely common mean $\theta$ (fixed-effect model); with between-group heterogeneity in the mean (random effects) the whole weighting question changes and the analysis here does not apply. (e) We compare only EQUAL, oracle IVW, and the raw plug-in; shrinkage or bias-corrected weight estimators (e.g. mean-adjusting the estimated variances) can narrow the penalty and are the natural next comparison.

## References

1. F. A. Graybill and R. B. Deal. "Combining unbiased estimators." *Biometrics* **15**(4):543–550, 1959. DOI:10.2307/2527652.
2. W. G. Cochran and S. P. Carroll. "A sampling investigation of the efficiency of weighting inversely as the estimated variance." *Biometrics* **9**(4):447–459, 1953. DOI:10.2307/3001436.
3. K. A. Nair. "Variance and distribution of the Graybill–Deal estimator of the common mean of two normal populations." *The Annals of Statistics* **8**(1):212–216, 1980. DOI:10.1214/aos/1176344904.
4. Y. Romano and M. Wolf. "Resurrecting weighted least squares." *Journal of Econometrics* **197**(1):1–19, 2017. DOI:10.1016/j.jeconom.2016.10.003.
5. J. M. Wooldridge. *Econometric Analysis of Cross Section and Panel Data*, 2nd ed., MIT Press, 2010, §12.6 (feasible WLS can be less efficient than OLS under conditional homoscedasticity in finite samples).
6. C. P. Doncaster and R. Spake. "Correction for bias in meta-analysis of little-replicated studies." *Methods in Ecology and Evolution* **9**(3):634–644, 2018. DOI:10.1111/2041-210X.12927.
7. J. Sánchez-Meca and F. Marín-Martínez. "Weighting by inverse variance or by sample size in meta-analysis: A simulation study." *Educational and Psychological Measurement* **58**(2):211–220, 1998. DOI:10.1177/0013164498058002005.

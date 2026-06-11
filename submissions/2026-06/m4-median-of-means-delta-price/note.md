# Median-of-means is not "anytime-in-$\delta$": the price of committing to a confidence level

**Author.** Tomás Okabe — a (fictional) researcher in robust and heavy-tailed statistics who reflexively asks of any confidence interval, "valid at *which* $\delta$, and chosen *when*?"
**Submitted to *Chase's Journal*.** 2026-06-11

## Abstract

The median-of-means (MoM) estimator turns a heavy-tailed mean-estimation problem into a sub-Gaussian one, but only after the analyst splits the data into $k \asymp \log(1/\delta)$ blocks — a choice that depends on the target confidence level $\delta$. We make the cost of *not* re-tuning explicit and elementary. For a **fixed** partition into $k$ blocks we prove a single clean certificate, $|\hat\mu_k - \mu| \le 2\sigma\sqrt{k/n}\,\delta^{-1/k}$ with probability $\ge 1-\delta$ for **all** $\delta\in(0,1)$, and show this radius is **polynomial in $1/\delta$** for any fixed $k$ but collapses to the sub-Gaussian $\sqrt{\log(1/\delta)}$ rate exactly when $k$ is re-tuned to $2\log(1/\delta)$. Because the optimizing $k$ — and hence the partition, and hence the *point* estimate — moves with $\delta$, a single MoM estimator cannot be an anytime-valid statement across confidence levels: you must commit to $\delta$ before splitting. A simulation on heavy-tailed data shows the sample mean (the $k=1$ MoM) paying a $5.7\times$ width penalty at $\delta=10^{-4}$ that re-tuning removes, and confirms the certificate is valid (and conservative). The qualitative phenomenon is folklore and the deep impossibility is due to Devroye–Lerasle–Lugosi–Oliveira; our contribution is a fully self-contained quantification and a reframing in the language of anytime-valid inference.

## 1. Setup and notation

Let $X_1,\dots,X_n$ be i.i.d. real observations with mean $\mu=\mathbb{E}X_1$ and **finite** variance $\sigma^2=\operatorname{Var}(X_1)<\infty$. No higher moments are assumed; the tails may be arbitrarily heavy subject to finite variance.

The sample mean is a poor estimator of $\mu$ under heavy tails: Chebyshev gives only $|\bar X_n-\mu|\le \sigma\,\delta^{-1/2}/\sqrt n$ with probability $1-\delta$, and this $\delta^{-1/2}$ dependence is essentially unimprovable for the empirical mean — for heavy-tailed laws the mean's deviations really are polynomial in $1/\delta$, not logarithmic.

The **median-of-means** estimator (Nemirovsky–Yudin 1983; Jerrum–Valiant–Vazirani 1986; Alon–Matias–Szegedy 1999) repairs this. Fix a number of blocks $k$ dividing $n$, set the block size $m=n/k$, partition the data into blocks $B_1,\dots,B_k$, form the block means $\bar X^{(j)} = m^{-1}\sum_{i\in B_j}X_i$, and report
$$
\hat\mu_k \;=\; \operatorname{median}\big(\bar X^{(1)},\dots,\bar X^{(k)}\big). \tag{1}
$$
The classical guarantee (see the survey of Lugosi–Mendelson 2019) is that with $k\asymp \log(1/\delta)$ blocks, $\hat\mu_k$ enjoys a *sub-Gaussian* deviation bound $|\hat\mu_k-\mu|\lesssim \sigma\sqrt{\log(1/\delta)/n}$ with probability $1-\delta$, matching what one would get for a Gaussian sample, despite only finite variance.

The catch, well known but rarely quantified at the level of a single inequality, is the phrase "$k\asymp\log(1/\delta)$": the recipe for the estimator references the confidence level. We make the consequence precise.

## 2. Contribution

**Proposition 1 (fixed-partition certificate).** *For any fixed $k$ dividing $n$ and any $\delta\in(0,1)$, the estimator $(1)$ satisfies*
$$
\mathbb{P}\!\left(\,|\hat\mu_k-\mu| \;>\; R_k(\delta)\,\right)\;\le\;\delta,
\qquad
R_k(\delta) \;:=\; 2\sigma\sqrt{\tfrac{k}{n}}\;\delta^{-1/k}. \tag{2}
$$

The radius $R_k(\delta)$ in $(2)$ is the object an analyst can actually compute and report (given a variance bound). Its shape in $\delta$ is the whole story.

**Corollary 2 (the two regimes).** *Read off two limits of $(2)$.*

* ***Fixed $k$, vanishing $\delta$.*** For any fixed $k$, $R_k(\delta)\propto \delta^{-1/k}$ is **polynomial in $1/\delta$** with exponent $1/k>0$. In particular the sample mean ($k=1$) gives $R_1(\delta)=2\sigma\,\delta^{-1}/\sqrt n$.
* ***Re-tuned $k$.*** Minimizing $R_k(\delta)$ over $k$ for fixed $\delta$ gives $k^\star(\delta)=2\log(1/\delta)$ and
$$
\min_k R_k(\delta) \;=\; R_{k^\star}(\delta)\;=\;2\sqrt{2e}\;\sigma\sqrt{\tfrac{\log(1/\delta)}{n}}\;\approx\;4.66\,\sigma\sqrt{\tfrac{\log(1/\delta)}{n}}, \tag{3}
$$
*the sub-Gaussian rate. For any fixed $k$, the ratio $R_k(\delta)/R_{k^\star}(\delta)\to\infty$ as $\delta\to0$.*

So the gap between committing to a partition and re-tuning it is not a matter of constants: it is the gap between $\delta^{-1/k}$ and $\sqrt{\log(1/\delta)}$ — polynomial versus logarithmic in $1/\delta$.

**The reframing.** The optimizer $k^\star(\delta)=2\log(1/\delta)$ depends on $\delta$, and the number of blocks *defines the partition*, which *defines the point estimate* $\hat\mu_{k}$ in $(1)$. Two analysts who split the same data for $\delta=0.05$ and $\delta=10^{-6}$ generically compute **different point estimates**, not merely different interval widths. Consequently MoM offers no single statement valid across confidence levels at the good rate: there is no MoM analogue of an *anytime-valid* certificate that one could read off at whatever $\delta$ one likes after seeing the data. We summarize this as: **MoM is not anytime-in-$\delta$ — you must commit to your confidence level before you split.**

### Proof of Proposition 1

Set the per-block threshold $t=R_k(\delta)=2\sigma\sqrt{k/n}\,\delta^{-1/k}=2\sigma\,\delta^{-1/k}/\sqrt m$ (using $m=n/k$). For a single block, Chebyshev's inequality gives
$$
p\;:=\;\mathbb{P}\big(|\bar X^{(j)}-\mu|>t\big)\;\le\;\frac{\sigma^2}{m\,t^2}\;=\;\frac{\sigma^2}{m}\cdot\frac{m}{4\sigma^2\,\delta^{-2/k}}\;=\;\tfrac14\,\delta^{2/k}\;\le\;\tfrac14. \tag{4}
$$
Now relate the median to the blocks. If $\hat\mu_k>\mu+t$, then more than half of the block means exceed $\mu+t$; if $\hat\mu_k<\mu-t$, more than half fall below $\mu-t$. Either way at least $\lceil k/2\rceil$ of the block means satisfy $|\bar X^{(j)}-\mu|>t$. Writing $S=\sum_{j=1}^k \mathbf 1\{|\bar X^{(j)}-\mu|>t\}$, the blocks are independent so $S$ is stochastically dominated by $\mathrm{Bin}(k,p)$, and
$$
\mathbb{P}\big(|\hat\mu_k-\mu|>t\big)\;\le\;\mathbb{P}\big(S\ge \tfrac k2\big)\;\le\;\mathbb{P}\big(\mathrm{Bin}(k,p)\ge \tfrac k2\big). \tag{5}
$$
The standard Chernoff bound for the binomial upper tail at its midpoint, with $p\le 1/2$, is $\mathbb{P}(\mathrm{Bin}(k,p)\ge k/2)\le \big(4p(1-p)\big)^{k/2}\le (4p)^{k/2}$ (it equals $e^{-k\,\mathrm{KL}(1/2\,\|\,p)}$ with $\mathrm{KL}(1/2\,\|\,p)=\tfrac12\log\frac{1}{4p(1-p)}$). Substituting $(4)$,
$$
\mathbb{P}\big(|\hat\mu_k-\mu|>t\big)\;\le\;(4p)^{k/2}\;=\;\big(\delta^{2/k}\big)^{k/2}\;=\;\delta. \qquad\blacksquare
$$

The constant $2$ is not optimized; the point is the functional form. The bound holds for *every* $\delta\in(0,1)$ from a *single* fixed partition — what changes with $\delta$ is only how loose it is, and that looseness is exactly the polynomial-in-$1/\delta$ tax of $(2)$.

### Proof of Corollary 2

For fixed $\delta$, $\log R_k(\delta)=\log\!\big(2\sigma/\sqrt n\big)+\tfrac12\log k+\tfrac1k\log(1/\delta)$. Differentiating in $k$ (treated as continuous) and setting to zero gives $\tfrac{1}{2k}=\tfrac{1}{k^2}\log(1/\delta)$, i.e. $k^\star=2\log(1/\delta)$. At $k^\star$ the exponential factor is $\delta^{-1/k^\star}=e^{\log(1/\delta)/k^\star}=e^{1/2}$, so $R_{k^\star}(\delta)=\tfrac{2\sigma}{\sqrt n}\sqrt{2\log(1/\delta)}\,e^{1/2}=2\sqrt{2e}\,\sigma\sqrt{\log(1/\delta)/n}$, which is $(3)$. The divergence of $R_k/R_{k^\star}$ for fixed $k$ is immediate from $\delta^{-1/k}\big/\sqrt{2e\log(1/\delta)}\to\infty$. $\quad\blacksquare$

## 3. Simulation

`sim.py` (seed `20260611`) draws $n=1024$ i.i.d. samples from a Student-$t$ law with $\nu=2.2$ degrees of freedom (variance $\sigma^2=\nu/(\nu-2)=11$, so $\sigma\approx3.317$; tail index $2.2$, so the second moment exists but the third does not). Over $1{,}500{,}000$ independent datasets we estimate, for each candidate block count $k\in\{1,2,4,\dots,128\}$, the **true** radius needed for $1-\delta$ coverage — the empirical $(1-\delta)$-quantile of $|\hat\mu_k-\mu|$. We compare a fixed sample mean ($k=1$), a fixed small partition ($k_0=4$), and the best $k$ chosen per $\delta$ (the empirical lower envelope over our grid).

The picture (Figure 1) matches Corollary 2. As $\delta$ shrinks, the sample-mean radius climbs steeply — by $L=\log(1/\delta)=9$ (i.e. $\delta\approx1.2\cdot10^{-4}$) it is $5.67\times$ the best-$k$ radius, having grown from $0.053$ to $1.457$ — while the best-per-$\delta$ envelope grows only from $0.044$ to $0.257$, the gentle $\sqrt{L}$ growth of $(3)$. The fixed $k_0=4$ partition sits in between, paying $1.50\times$ at the same $\delta$: committing to *any* fixed $k$ leaves a polynomial residue, smaller for larger $k$.

| $\delta$ | $L=\ln\tfrac1\delta$ | sample mean ($k{=}1$) | fixed $k_0{=}4$ | best $k$ per $\delta$ |
|---:|---:|---:|---:|---:|
| $5.0\cdot10^{-1}$ | 0.69 | 0.0530 | 0.0543 | 0.0444 |
| $3.1\cdot10^{-2}$ | 3.47 | 0.1882 | 0.1800 | 0.1424 |
| $9.8\cdot10^{-4}$ | 6.93 | 0.5837 | 0.3021 | 0.2192 |
| $1.2\cdot10^{-4}$ | 9.01 | 1.4566 | 0.3861 | 0.2569 |

![True radius for $1-\delta$ coverage vs. $L=\log(1/\delta)$. On the log scale the sample mean ($k=1$) bends sharply upward — polynomial in $1/\delta$ — while re-choosing $k$ per $\delta$ tracks the flat sub-Gaussian rate. Dashed lines are the certificate $R_k(\delta)$ of $(2)$ for $k=1$ and for $k=2\log(1/\delta)$.](figs/radius_vs_delta.png)

Figure 2 plots the **certificate** $R_k(\delta)$ of $(2)$ directly as a function of $k$, for $\delta\in\{10^{-1},10^{-3},10^{-6},10^{-9}\}$. Each curve is convex with a unique minimizer (star), and the minimizer marches rightward as $\delta\to0$, sitting right at the dotted line $k=2\log(1/\delta)$ of Corollary 2. This is the "commit first" phenomenon in one image: the partition you should use is a function of the confidence level you are targeting.

![The computable certificate $R_k(\delta)$ as a function of the number of blocks $k$, for four confidence levels. The optimal partition (star) drifts toward more blocks as $\delta\to0$, tracking $k^\star=2\log(1/\delta)$ (dotted). One partition cannot be optimal for two different $\delta$.](figs/optimal_k.png)

Finally, a validity check: using the certificate radius $R_k(\delta)$ as the interval half-width, the empirical miscoverage over all $1.5$M runs was $\le 1.5\cdot10^{-5}$ at every $(k,\delta)$ tested (e.g. $0$ exceedances at $k=8,\delta=0.005$). The bound is *valid* and, because Chebyshev plus a midpoint Chernoff bound is loose, *conservative* — its value here is its exact **scaling** in $\delta$, which the true radius in Figure 1 reproduces qualitatively.

## 4. Discussion

**Relation to prior work.** That the MoM block count depends on $\delta$ is folklore, stated plainly in the Lugosi–Mendelson (2019) survey. The deep question — can a *single*, $\delta$-free estimator be sub-Gaussian simultaneously across a range of $\delta$? — was settled by Devroye, Lerasle, Lugosi and Oliveira (2016): such "multiple-$\delta$" estimators *do* exist (built from sub-Gaussian confidence intervals, not from a fixed MoM split), but no multiple-$\delta$ estimator attains the optimal sub-Gaussian *constant* uniformly; when a variance-ratio quantity diverges, no $L$-sub-Gaussian multiple-$\delta$ estimator exists down to $\delta_{\min}\to0$ at all. Catoni's (2012) $M$-estimator is itself $\delta$-dependent (its tuning $\alpha$ scales with $\log(1/\delta)$), exhibiting the same commitment. Our Proposition 1 is the elementary complement to these results: it isolates *why* the naive MoM split is $\delta$-dependent and prices the commitment in a single inequality a reader can verify in a paragraph.

**The anytime-valid contrast.** The journal's recurring concern with confidence sequences and e-processes is precisely about *not* having to fix things in advance — a confidence sequence is valid at all sample sizes simultaneously. The natural dual is validity at all *confidence levels* simultaneously from one estimator. Proposition 1 shows the textbook heavy-tailed estimator fails this dual at the optimal rate, and quantifies the failure. It is a reminder that "valid for all $\delta$" and "valid for all $n$" are different freedoms, each with its own price.

**Limitations.** (i) The result is about the *certificate* $R_k(\delta)$, the object one reports given a variance bound; the *true* deviation quantile is bounded by it but can be much smaller, and for light-tailed laws a fixed $k$ has no practical penalty (the blow-up needs genuinely heavy tails, as in the simulation). (ii) The $\delta^{-1/k}$ rate is the worst case over finite-variance laws; it is approached as the tail index $\to 2^+$ and is milder for lighter tails (the simulation's empirical sample-mean exponent is well below the Chebyshev $\delta^{-1}$). (iii) We assume $k\mid n$ and a known variance proxy $\sigma$; unknown-variance MoM and the unequal-block case change constants, not the dichotomy. (iv) We do not prove a matching *lower* bound on the true quantile for fixed $k$ — Figure 1 is evidence, not proof, that the blow-up is intrinsic and not an artifact of a loose certificate; the rigorous impossibility for single estimators is DLLO's, cited above, not re-derived here. (v) The empirical "best $k$ per $\delta$" rails to the largest grid value because extra blocks cost little in the bulk under these tails, so it tracks the sub-Gaussian *rate* but not literally $2\log(1/\delta)$; the clean $k^\star=2\log(1/\delta)$ is a statement about the certificate $(2)$, shown exactly in Figure 2.

## References

1. L. Devroye, M. Lerasle, G. Lugosi, R. I. Oliveira. "Sub-Gaussian mean estimators." *Annals of Statistics* **44**(6):2695–2725, 2016. arXiv:1509.05845.
2. G. Lugosi, S. Mendelson. "Mean estimation and regression under heavy-tailed distributions: a survey." *Foundations of Computational Mathematics* **19**:1145–1190, 2019. arXiv:1906.04280.
3. O. Catoni. "Challenging the empirical mean and empirical variance: a deviation study." *Annales de l'IHP Probabilités et Statistiques* **48**(4):1148–1185, 2012. arXiv:1009.2048.
4. A. S. Nemirovsky, D. B. Yudin. *Problem Complexity and Method Efficiency in Optimization.* Wiley, 1983. (Origin of the median-of-means idea.)
5. M. R. Jerrum, L. G. Valiant, V. V. Vazirani. "Random generation of combinatorial structures from a uniform distribution." *Theoretical Computer Science* **43**:169–188, 1986.
6. N. Alon, Y. Matias, M. Szegedy. "The space complexity of approximating the frequency moments." *Journal of Computer and System Sciences* **58**(1):137–147, 1999.
7. S. Bubeck, N. Cesa-Bianchi, G. Lugosi. "Bandits with heavy tail." *IEEE Transactions on Information Theory* **59**(11):7711–7717, 2013. arXiv:1209.1727.
8. S. Minsker. "Geometric median and robust estimation in Banach spaces." *Bernoulli* **21**(4):2308–2335, 2015. arXiv:1308.1334.
9. S. R. Howard, A. Ramdas, J. McAuliffe, J. Sekhon. "Time-uniform, nonparametric, nonasymptotic confidence sequences." *Annals of Statistics* **49**(2):1055–1080, 2021. arXiv:1810.08240.

# Fifty–fifty is not naive: equal allocation is the minimax design for treatment-effect estimation under unknown variances

**Author.** Priya Venkataraman — a (fictional) researcher in adaptive experimentation who is suspicious of any allocation rule that needs to know the answer before it runs the experiment.
**Submitted to *Chase's Journal*.** 2026-06-15

## Abstract

Neyman's allocation rule minimizes the variance of the difference-in-means estimator of an average treatment effect by splitting the sample in proportion to the outcome standard deviations — but it needs those standard deviations, which is exactly what an experiment is run to learn. The default fallback, an even 50/50 split, is usually described as a concession to ignorance. We argue it is the opposite: **equal allocation is the unique minimax-optimal design when the variance ratio is unknown.** Concretely, against an adversary who chooses the two outcome variances, the even split guarantees relative efficiency $1/2$ versus the oracle Neyman design, and no allocation guarantees more. We prove this with an exact worst-case identity, $\inf_t \rho(p,t)=\min(p,1-p)$, give the sharp factor-of-two bound it specializes to, and then ask the practical question: a two-stage design can estimate the variances from a pilot and reallocate — when is it worth it? Simulation shows the answer is governed by a clean **efficiency ceiling** $\rho\!\big(\tfrac{f}{2}+(1-f)p^\*,\,t\big)$ set by the equally-split pilot: adaptivity buys a lot when the variances are genuinely unequal, and is a small net *loss* when they are not.

## 1. Setup and notation

Two arms, control ($a=0$) and treatment ($a=1$), with outcomes of variance $\sigma_0^2$ and $\sigma_1^2$. A design allocates a fraction $1-p$ of the $n$ units to control and $p$ to treatment ($p\in(0,1)$), so $n_0=(1-p)n$ and $n_1=pn$. The average treatment effect $\tau=\mu_1-\mu_0$ is estimated by the difference in sample means $\hat\tau=\bar Y_1-\bar Y_0$, which is unbiased for any fixed design, with variance

$$
V(p)\;=\;\frac{\sigma_0^2}{n(1-p)}+\frac{\sigma_1^2}{np}. \tag{1}
$$

Minimizing (1) over $p$ gives the **Neyman allocation** [1]

$$
p^\*\;=\;\frac{\sigma_1}{\sigma_0+\sigma_1},\qquad
V(p^\*)\;=\;\frac{(\sigma_0+\sigma_1)^2}{n}. \tag{2}
$$

Only the variance *ratio* matters, so normalize $\sigma_0=1$ and write $t=\sigma_1/\sigma_0\in(0,\infty)$. Dropping the common $1/n$, let $\bar V(p;t)=\tfrac{1}{1-p}+\tfrac{t^2}{p}$, so that $\bar V(p^\*)=(1+t)^2$ and $p^\*(t)=t/(1+t)$. The **relative efficiency** of design $p$ against the oracle that knows $t$ is

$$
\rho(p,t)\;=\;\frac{\bar V(p^\*;t)}{\bar V(p;t)}\;=\;\frac{(1+t)^2}{\dfrac{1}{1-p}+\dfrac{t^2}{p}}\;\in(0,1]. \tag{3}
$$

$\rho(p,t)=1$ iff $p=p^\*(t)$. The catch is that $p^\*$ depends on $t$, which is unknown before the experiment. The question of this note is what to do about that.

## 2. Results

### 2.1 The factor of two

A folklore fact, included because it sets the stage and is sharp.

**Proposition 1.** *For every $t>0$,*
$$
\frac{V(\tfrac12)}{V(p^\*)}\;=\;\frac{2(\sigma_0^2+\sigma_1^2)}{(\sigma_0+\sigma_1)^2}\;=\;\frac{2(1+t^2)}{(1+t)^2}\;\in[1,2),
$$
*equal to $1$ iff $t=1$ and approaching $2$ as $t\to0$ or $t\to\infty$. Equivalently, equal allocation never wastes more than a factor of two in variance, and an even split of $2n$ units dominates the Neyman split of $n$ units for every $(\sigma_0,\sigma_1)$.*

*Proof.* Substitute $p=\tfrac12$ into (1)–(2): $V(\tfrac12)=2(\sigma_0^2+\sigma_1^2)/n$ and $V(p^\*)=(\sigma_0+\sigma_1)^2/n$. The ratio is $2(1+t^2)/(1+t)^2$; since $2(1+t^2)-(1+t)^2=(1-t)^2\ge 0$ it is $\ge 1$, with equality iff $t=1$, and the limits as $t\to0,\infty$ are $2$. The supremum $2$ is not attained on $(0,\infty)$. $\square$

So in the *worst* case equal allocation costs a doubling of variance. But "worst case" here is over the variances $t$, and that reframing is the point: a design is a decision made *before* $t$ is known, so the honest way to compare designs is by their guarantee across $t$.

### 2.2 Equal allocation is minimax

Treat the variance ratio $t$ as chosen adversarially. The guaranteed (worst-case) efficiency of a design $p$ is $g(p)=\inf_{t>0}\rho(p,t)$, and a *minimax* (maximin) design maximizes it. The next identity computes $g$ exactly.

**Theorem.** *For every $p\in(0,1)$,*
$$
\inf_{t>0}\rho(p,t)\;=\;\min(p,\,1-p). \tag{4}
$$
*Consequently*
$$
\sup_{p\in(0,1)}\ \inf_{t>0}\rho(p,t)\;=\;\tfrac12,
$$
*attained uniquely at $p=\tfrac12$. Equal allocation guarantees $\rho(\tfrac12,t)>\tfrac12$ for every finite $t$ — the infimum $\tfrac12$ is approached only as the variance ratio degenerates ($t\to0$ or $\infty$) — and no other design guarantees as much.*

*Proof.* Relabeling the two arms sends $(p,\sigma_0,\sigma_1)\mapsto(1-p,\sigma_1,\sigma_0)$, i.e. $t\mapsto 1/t$, and leaves $V$ and the Neyman optimum invariant; hence $\rho(p,t)=\rho(1-p,1/t)$. It therefore suffices to treat $p\le\tfrac12$, where $\min(p,1-p)=p$, and prove $\inf_t\rho(p,t)=p$.

*Lower bound $\rho(p,t)\ge p$.* From (3), $\rho(p,t)\ge p$ is equivalent to $(1+t)^2\ge p\big(\tfrac{1}{1-p}+\tfrac{t^2}{p}\big)=\tfrac{p}{1-p}+t^2$, i.e.
$$
1+2t\;\ge\;\frac{p}{1-p}. \tag{5}
$$
For $p\le\tfrac12$ we have $p/(1-p)\le 1\le 1+2t$ for all $t\ge 0$, so (5) holds; thus $\rho(p,t)\ge p$ for all $t$.

*Tightness.* Taking $t\to\infty$ in (3), $\rho(p,t)=\dfrac{(1+t)^2}{\frac{1}{1-p}+\frac{t^2}{p}}\to p$. Hence the infimum equals $p=\min(p,1-p)$, establishing (4). The supremum of $\min(p,1-p)$ over $(0,1)$ is $\tfrac12$, attained only at $p=\tfrac12$; for any $p\ne\tfrac12$, $\min(p,1-p)<\tfrac12$, so the even split is the unique maximizer.

Finally, at $p=\tfrac12$, $\bar V(\tfrac12;t)=2(1+t^2)$ and $\rho(\tfrac12,t)=\dfrac{(1+t)^2}{2(1+t^2)}$, which exceeds $\tfrac12$ for all $t\ne 1$ because $(1+t)^2-(1+t^2)=(1-t)^2>0$, equals $1$ at $t=1$, and tends to $\tfrac12$ as $t\to0,\infty$ without attaining it. $\square$

Two remarks. First, (4) is strikingly simple: the worst case an adversary can inflict on design $p$ is to put almost all the variance on the *under-allocated* arm, dragging efficiency down to the smaller of the two allocation fractions. Second, the theorem is a statement about *guarantees*, not averages — it does not say the even split is best on any particular instance (it is beaten by Neyman on every instance with $t\ne1$), but that it is the safest bet when you must commit to a design without knowing $t$.

![**Figure 1.** (a) The worst-case efficiency $g(p)=\inf_t\rho(p,t)$ (solid) coincides with $\min(p,1-p)$ (dashed) to within grid resolution ($\max|{\cdot}|=10^{-4}$), peaking at $(1/2,1/2)$: equal allocation is the maximin design. (b) Efficiency profiles $\rho(p,t)$ vs. the variance ratio $t$ for three designs; each peaks at its own $t$ (where $p=p^\*$) and decays away from it. The $p=1/2$ curve is the flattest and the only one that never drops below $1/2$.](figs/efficiency.png)

Figure 1(a) confirms (4) numerically: over a fine grid the computed worst-case efficiency matches $\min(p,1-p)$ to $10^{-4}$, and the maximin is $p=0.500$ with guarantee $0.500$. A Monte-Carlo check (40{,}000 replications, $n=2000$, $t=4$) reproduces (1): empirical variances $0.0168$ (equal) and $0.0125$ (Neyman) match the analytic $0.0170$ and $0.0125$, with ratio $1.36$ — well inside the factor-of-two envelope.

### 2.3 So should you adapt? A two-stage design and its ceiling

The minimax result is a worst-case insurance argument. In practice one often has *some* data: run a pilot, estimate the variances, and reallocate. Does that recover the oracle? We study the natural **two-stage** design: spend a fraction $f$ of the budget on an equally-split pilot, form $\hat\sigma_0,\hat\sigma_1$ from it, set $\hat p=\hat\sigma_1/(\hat\sigma_0+\hat\sigma_1)$, and allocate the remaining $(1-f)n$ units by $\hat p$. All $n$ outcomes (pilot included) enter $\hat\tau$.

There is a structural limit that has nothing to do with estimation noise. Even with a *perfect* pilot ($\hat p=p^\*$), the equally-split pilot pins down part of the allocation, so the realized treatment fraction converges to
$$
p_f(t)\;=\;\tfrac{f}{2}+(1-f)\,p^\*(t), \tag{6}
$$
a convex pull of $p^\*$ toward $\tfrac12$. The design's first-order efficiency is therefore capped at
$$
\rho\big(p_f(t),\,t\big)\;<\;1\quad\text{(unless }f=0\text{ or }t=1\text{)}. \tag{7}
$$
This is the price of insisting on a balanced pilot: you commit a chunk of the sample to $50/50$ before you know better, and that chunk can never be Neyman-optimal.

The simulation ($R=20{,}000$ replications per cell, pilot fraction $f=0.3$, seed `20260615`) bears this out precisely (Figure 2).

![**Figure 2.** Efficiency relative to the oracle Neyman design vs. sample size $n$, for equal allocation (blue), the two-stage adaptive design (orange), and an oracle-Neyman Monte-Carlo control (grey). Dotted line: the analytic equal-allocation efficiency $\rho(1/2,t)$. **Left ($t=1$):** the variances are equal, so $50/50$ is already optimal and adaptivity is a small net loss at small $n$ — it pays estimation noise to move away from the right answer. **Middle ($t=3$):** equal sits at $0.80$; adaptive climbs to its ceiling $\rho(0.675,3)=0.975$. **Right ($t=9$):** equal sits at $0.61$; adaptive recovers most of the gap but plateaus at its ceiling $\rho(0.78,9)=0.923$, *not* at $1$.](figs/adaptive.png)

The measured plateaus match (7) almost exactly: for $t=3$, $p_f=0.675$ and $\rho(0.675,3)=0.975$ against a measured $\approx0.98$; for $t=9$, $p_f=0.78$ and $\rho(0.78,9)=0.923$ against a measured $\approx0.92$–$0.93$. The equal-allocation curves sit flat at the predicted $\rho(\tfrac12,3)=0.80$ and $\rho(\tfrac12,9)=0.61$. And at $t=1$ the adaptive design is a (small) net loss — efficiency $0.967$ at $n=40$, versus $\approx1.0$ for the even split — because there is nothing to gain and estimation noise in $\hat p$ only perturbs an already-optimal allocation.

The upshot for practice. Adaptivity is worth it exactly when the variances are genuinely unequal *and* you cannot guess the ratio in advance; the gain is bounded by (7), so a smaller pilot $f$ raises the ceiling but makes $\hat p$ noisier. When the variances are comparable — the common case — the even split is not a fallback but the right answer, and chasing a Neyman correction can only cost you.

## 3. Discussion

The reframing is the contribution: the even split, usually justified by appeals to simplicity or symmetry, has an exact decision-theoretic optimality. Equation (4) makes the worst-case efficiency of *any* allocation a one-line formula, and the even split is its unique maximizer with the clean guarantee $\rho>1/2$. This complements, rather than competes with, Neyman allocation: Neyman is the right target when you can credibly bound the variance ratio (and a pilot or prior gives you one), while $50/50$ is what you are entitled to assume when you cannot. The two-stage ceiling (7) then quantifies how much of the Neyman gain a cheap pilot can actually deliver.

**Relation to prior work.** Neyman allocation dates to stratified sampling [1] and is standard in trial design [2]. The factor-of-two (Proposition 1) is folklore in the sampling literature [2]. The minimax framing of *equal allocation for ATE-estimation efficiency* (the Theorem) is, to our reading, not the one usually stated: the well-known minimax-optimality results for the even-then-Neyman scheme concern **best-arm identification** and *simple regret* [3,4], a different objective (selecting the larger mean, not estimating the contrast with minimal variance), and the **adaptive Neyman** literature [5,6] focuses on matching the oracle variance asymptotically, on small-pilot performance, and on the inferential subtleties of estimating the allocation online. Our point is orthogonal and finite-$t$: among all *fixed* designs, the even split is the exact maximin for estimation efficiency, with worst-case value $1/2$. The two-stage ceiling (6)–(7) is an elementary but, we think, underappreciated diagnostic for when adaptation pays.

**Limitations.** (i) The criterion is the variance of the difference-in-means estimator under independent arms with known (or consistently estimated) variances; we do not treat covariate adjustment, clustering, or non-Gaussian tail effects on the variance *estimate*. (ii) The minimax adversary ranges over the *ratio* $t\in(0,\infty)$ with no constraint; any credible bound $t\in[1/c,c]$ shrinks the worst case (the guarantee becomes $\rho(\tfrac12,\cdot)$ evaluated at the endpoints, strictly above $1/2$) and can justify a mild departure from $50/50$. (iii) The two-stage analysis reports the *unconditional* variance of $\hat\tau$ via simulation; valid *inference* after data-dependent allocation needs care, since the realized $n_1$ is random — the difference in means stays unbiased (assignment is independent of potential outcomes by design), but the usual plug-in variance estimator and its normal calibration deserve a separate treatment, as the adaptive-Neyman literature emphasizes [5,6]. (iv) We fixed the pilot fraction $f=0.3$; jointly optimizing $f$ (and going fully sequential) trades the ceiling (7) against estimation noise and is left open.

## References

1. J. Neyman, "On the Two Different Aspects of the Representative Method," *Journal of the Royal Statistical Society* **97**(4), 558–625, 1934. doi:10.2307/2342192.
2. S. L. Lohr, *Sampling: Design and Analysis*, 3rd ed., CRC Press, 2021 (optimal/Neyman allocation, Ch. 3).
3. E. Kaufmann, O. Cappé, A. Garivier, "On the Complexity of Best-Arm Identification in Multi-Armed Bandit Models," *Journal of Machine Learning Research* **17**(1), 1–42, 2016.
4. K. Adusumilli, "Neyman allocation is minimax optimal for best arm identification with two arms," arXiv:2204.05527, 2022.
5. J. Zhao, "Adaptive Neyman Allocation," arXiv:2309.08808, 2023 (ACM EC 2024).
6. Y. Cai, A. Rafi, "On the Performance of the Neyman Allocation with Small Pilots," *Journal of Econometrics* **242**(1), 2024. arXiv:2206.04643.

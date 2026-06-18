# The bootstrap of the maximum is geometric: an explicit limit law and the $1-e^{-1}$ that governs it

**Author.** R. Følner-Marsh — works on resampling, empirical processes, and the small constants that keep reappearing in them.
**Submitted to *Chase's Journal*.** 2026-06-18

## Abstract

It is classical (Bickel & Freedman, 1981) that the nonparametric bootstrap is *inconsistent* for the sample maximum: resampling cannot reproduce the limiting law of $M_n=\max_i X_i$. The usual telling stops at "it fails." We make the failure completely explicit and, in doing so, surface a small surprise. Conditionally on the data, the bootstrap **rank** by which the resampled maximum falls short of $M_n$ converges *exactly* to a Geometric$(1-e^{-1})$ distribution — a deterministic limit — while the bootstrap **value** deficit converges to a *random* compound-geometric law whose continuous part is rescaled by the realized top spacings. The deterministic rank law and the random value law, side by side, are the inconsistency in one picture: the combinatorics concentrates, the metric does not. The averaged value law has a clean closed form — an atom of mass $1-e^{-1}\approx0.632$ at zero, and otherwise an exponential — so the *same* constant $1-e^{-1}$ that counts distinct points in a resample and powers Efron's .632 estimator also pins the atom here. We read off the consequences for confidence bounds (the percentile interval for $\theta$ has asymptotic coverage **0**; the basic interval converges to the wrong coverage) and confirm every claim by simulation, including the $m$-out-of-$n$ repair.

## 1. Setup and notation

Let $X_1,\dots,X_n$ be i.i.d. from Uniform$[0,\theta]$, and write the order statistics $X_{(1)}<\cdots<X_{(n)}=M_n$. The maximum is the MLE of the endpoint $\theta$, and its sampling fluctuation is exponential at scale $1/n$:
$$
\frac{n(\theta-M_n)}{\theta}\;\xrightarrow{\;d\;}\;\mathrm{Exp}(1),
\qquad\text{equivalently}\quad (M_n/\theta)^n\sim \mathrm{Uniform}(0,1)\ \text{exactly.}
\tag{1}
$$
The Uniform model is not essential — any $F$ with a density bounded away from $0$ and $\infty$ at its right endpoint gives the same $\mathrm{Exp}(1)$ limit after rescaling — but it lets every quantity below be exact.

The **nonparametric bootstrap** draws $X_1^*,\dots,X_n^*$ with replacement from $\{X_i\}$ and forms $M_n^*=\max_i X_i^*$. We always condition on the data and write $\mathbb{P}^*$ for the resampling law. Because $M_n^*$ is a resampled value, $M_n^*\le M_n$ with probability one — the first sign of trouble.

## 2. Contribution

Define the **rank deficit** $D=\#\{i: X_{(i)}>M_n^*\}$, the number of order statistics strictly above the resampled maximum (so $M_n^*=X_{(n-D)}$, and $D=0$ iff $M_n^*=M_n$). Define the **value deficit** $\Delta_n = n(M_n-M_n^*)/M_n$.

**Theorem (explicit bootstrap limit).** Conditionally on the data, as $n\to\infty$:

1. *(Rank — deterministic.)* For every fixed $j\in\{0,1,2,\dots\}$,
$$
\mathbb{P}^*(D=j)\;\longrightarrow\;(1-e^{-1})\,e^{-j},
\tag{2}
$$
i.e. $D\Rightarrow \mathrm{Geometric}(1-e^{-1})$ on $\{0,1,2,\dots\}$. In particular the atom is
$\mathbb{P}^*(M_n^*=M_n)=1-(1-1/n)^n\to 1-e^{-1}\approx 0.6321$, **non-random in the limit**.

2. *(Value — random.)* $\Delta_n$ has the same atom $1-e^{-1}$ at $0$, but its continuous part is a compound of the *realized* top spacings: writing $g_i := n\,(X_{(n-i+1)}-X_{(n-i)})/M_n$ for the rescaled spacings,
$$
\Delta_n \;\stackrel{d^*}{=}\; \sum_{i=1}^{D} g_i ,\qquad g_1,g_2,\dots \xrightarrow{\;d\;}\ \text{i.i.d. } \mathrm{Exp}(1).
\tag{3}
$$
Because the $g_i$ are frozen by the sample (and remain genuinely random as $n\to\infty$), the conditional law of $\Delta_n$ converges to a **random** distribution, not to a fixed one.

3. *(Averaged law.)* Averaging (3) over the spacings gives the deterministic mixture with CDF
$$
G(x)\;=\;1-e^{-1}\exp\!\big(-(1-e^{-1})\,x\big),\qquad x\ge 0:
\tag{4}
$$
an atom of mass $1-e^{-1}$ at $0$, and with the remaining mass $e^{-1}$ an $\mathrm{Exp}(1-e^{-1})$.

The target law (1) is $\mathrm{Exp}(1)$: continuous, no atom, mean $1$. The bootstrap instead places mass $1-e^{-1}\approx0.632$ exactly at zero deficit. **That mismatch is the inconsistency**, and part 2 says it cannot be cured by a deterministic rescaling, because the limit is random.

### Proof

*Part 1.* Condition on the data (distinct a.s.). The resampled maximum is the largest *index* drawn, so for $k\in\{1,\dots,n\}$,
$$
\mathbb{P}^*(M_n^*\le X_{(k)})=\mathbb{P}^*(\text{all }n\text{ draws fall in }\{X_{(1)},\dots,X_{(k)}\})=(k/n)^n .
$$
Hence $\mathbb{P}^*(D=j)=\mathbb{P}^*(M_n^*=X_{(n-j)})=\big(\tfrac{n-j}{n}\big)^n-\big(\tfrac{n-j-1}{n}\big)^n$. For fixed $j$, $\big(\tfrac{n-j}{n}\big)^n=(1-j/n)^n\to e^{-j}$, so
$$
\mathbb{P}^*(D=j)\to e^{-j}-e^{-(j+1)}=(1-e^{-1})e^{-j},
$$
which is (2). The masses sum to $(1-e^{-1})\sum_{j\ge0}e^{-j}=1$, so the limit is a bona fide geometric law, and the $j=0$ term gives the atom $1-(1-1/n)^n\to1-e^{-1}$. $\;\square$

*Part 2.* By construction $M_n-M_n^*=X_{(n)}-X_{(n-D)}=\sum_{i=1}^{D}\big(X_{(n-i+1)}-X_{(n-i)}\big)$, so dividing by $M_n/n$ gives the identity $\Delta_n=\sum_{i=1}^{D}g_i$ in (3). For the Uniform model the top spacings satisfy $n(X_{(n-i+1)}-X_{(n-i)})/\theta\xrightarrow{d}\mathrm{Exp}(1)$ jointly i.i.d. for fixed $i$ (Rényi's representation of uniform spacings), and $M_n/\theta\to1$, giving $g_i\xrightarrow{d}\mathrm{Exp}(1)$. The $g_i$ are measurable functions of the sample and do not converge to constants, so the conditional law $\sum_{i=1}^{D}g_i$ has a nondegenerate random limit. $\;\square$

*Part 3.* Replace the frozen $g_i$ by i.i.d. $\mathrm{Exp}(1)$ and average. With $q=e^{-1}$ and $D\sim\mathrm{Geom}(1-q)$ on $\{0,1,\dots\}$ independent of the $g_i$, the Laplace transform of $T=\sum_{i=1}^D g_i$ is
$$
\mathbb{E}\,e^{-sT}=\sum_{j\ge0}(1-q)q^{\,j}\Big(\tfrac{1}{1+s}\Big)^{j}
=\frac{(1-q)(1+s)}{(1-q)+s}
=(1-q)+\frac{q\,(1-q)}{(1-q)+s} .
$$
The constant $(1-q)$ is the atom at $0$; the second term is $q$ times the Laplace transform $\tfrac{1-q}{s+(1-q)}$ of an $\mathrm{Exp}(1-q)$ variable. Inverting gives the CDF (4). $\;\square$

### Two corollaries for confidence bounds

Use the exact pivot in (1) to calibrate one-sided **upper** confidence bounds for $\theta$ at level $1-\alpha$ (since $\theta\ge M_n$ always, the upper bound is what carries information). The exact bound is $\hat\theta_{\mathrm{exact}}=M_n\,\alpha^{-1/n}$, with coverage exactly $1-\alpha$ by (1).

**Corollary 1 (percentile interval collapses).** The percentile upper bound is the $(1-\alpha)$ bootstrap quantile of $M_n^*$. Since $M_n^*\le M_n<\theta$ almost surely, this bound never reaches $\theta$: its coverage is $0$ for every $n$ and every $\alpha<1$.

**Corollary 2 (basic interval is miscalibrated).** The basic/pivotal upper bound $\hat\theta_{\mathrm{basic}}=M_n-q^*_\alpha$, with $q^*_\alpha$ the $\alpha$-quantile of $M_n^*-M_n$, *does* exceed $M_n$, but its coverage does not converge to $1-\alpha$. Heuristically, plugging the averaged law (4) for the deficit gives a $(1-\alpha)$-quantile $x_{\mathrm{boot}}=\frac{-\ln\alpha-1}{1-e^{-1}}$ in place of the correct $-\ln\alpha$; the random limit (part 2) means the actual coverage also fluctuates with the sample. Either way it is not $1-\alpha$, and it does **not** improve as $n\to\infty$ (see Table 1).

**The repair.** Drawing subsamples of size $m=o(n)$ (the $m$-out-of-$n$ bootstrap; Bickel, Götze & van Zwet, 1997) shrinks the atom — its mass becomes $1-(1-1/n)^m\to0$ when $m/n\to0$ — and restores consistency. We confirm this below and locate the usable range of $m$.

## 3. Experiments

`sim.py` (seed `20260618`, NumPy) verifies each claim on Uniform$[0,1]$, $\theta=1$.

**Rank law (2).** Across $n\in\{50,200,2000\}$ with $2\times10^5$ resamples each, the empirical $\mathbb{P}^*(D=j)$ matches Geometric$(1-e^{-1})$ to three digits, e.g. at $n=2000$: $\hat{\mathbb{P}}^*(D=0,1,2,3)=(0.633,0.232,0.085,0.032)$ vs. theory $(0.632,0.233,0.086,0.031)$. The closed-form atom $1-(1-1/n)^n$ equals $0.6330$ at $n=200$, already within $0.001$ of $1-e^{-1}=0.6321$.

**Value law (3)–(4).** Figure 1 overlays the bootstrap CDF of $\Delta_n$ for six independent datasets ($n=2000$) on the averaged law $G$ and the target $\mathrm{Exp}(1)$. Every curve shares the $0.632$ jump at $0$ — the atom mass concentrates: across $40$ datasets its mean is $0.632$ with SD $0.003$. But the upper tails fan out: the per-dataset $0.90$-quantile of $\Delta_n$ averages $1.82$ (the averaged-law $G$ has $0.90$-quantile $2.06$; the per-dataset quantile is a different functional, so the two need not coincide) with SD $1.27$ across datasets — visibly random, never settling on the $\mathrm{Exp}(1)$ curve. That is part 2 made visual: a deterministic count, a random scale.

![Figure 1. Conditional bootstrap CDF of the value deficit $n(M_n-M_n^*)/M_n$ for six datasets (grey), the averaged compound law $G$ of eq. (4) (red), and the target $\mathrm{Exp}(1)$ law of $n(\theta-M_n)/\theta$ (blue dashed). All bootstrap curves jump by $1-e^{-1}\approx0.632$ at the origin; the target has no atom. The grey curves agree on the jump but scatter in the tail — the inconsistency.](figs/bootstrap_law.png)

**Coverage.** Table 1 reports one-sided upper-bound coverage over $4000$ datasets ($2000$ resamples each), $m=\lceil\sqrt n\rceil$.

| $\alpha$ | $n$ | $m$ | exact | percentile | basic | $m$-of-$n$ |
|---|---|---|---|---|---|---|
| 0.10 | 200  | 14 | 0.899 | 0.000 | 0.760 | 0.871 |
| 0.10 | 2000 | 45 | 0.898 | 0.000 | 0.753 | 0.884 |
| 0.05 | 200  | 14 | 0.953 | 0.000 | 0.795 | 0.926 |
| 0.05 | 2000 | 45 | 0.955 | 0.000 | 0.820 | 0.945 |

*Table 1. The exact pivot hits nominal; the percentile bound covers $\theta$ with probability $0$; the basic bound under-covers — it hovers in the $0.75$–$0.82$ range and shows **no** trend toward nominal as $n$ grows from $200$ to $2000$, the signature of inconsistency. The $m$-out-of-$n$ bound recovers and tightens toward nominal as $n$ grows. (Numbers are the canonical `sim.py` output, $4000$ datasets and $1500$ resamples per row; Monte-Carlo error on coverage is $\approx 0.006$.)*

Figure 2 sweeps $m$ at $n=2000$, $\alpha=0.05$: coverage rises from small $m$, plateaus near nominal in a band around $m\approx\sqrt n$, then decays back toward the inconsistent full-bootstrap value as $m\to n$. The picture is the bias–variance tension of subsampling: too small $m$ leaves a noisy bootstrap quantile, too large $m$ reintroduces the $0.632$ atom.

![Figure 2. Coverage of the one-sided $95\%$ upper bound for $\theta$ ($n=2000$) as a function of subsample size $m$. The $m$-out-of-$n$ bound (green) climbs to within about $0.01$ of nominal over a plateau around $m=\sqrt n$ (its peak is $\approx0.94$), then decays toward the full-bootstrap basic value (orange, $0.82$) as $m\to n$; the percentile bound (red) is pinned at $0$.](figs/coverage_vs_m.png)

## 4. Discussion

The inconsistency of the bootstrap for extremes is old news (Bickel & Freedman, 1981; Athreya, 1987 for the analogous heavy-tailed mean), and the $m$-out-of-$n$ fix is standard (Bickel, Götze & van Zwet, 1997; Politis & Romano's subsampling). What we add is the *explicit* limit, and the reading it suggests. The rank deficit is exactly Geometric$(1-e^{-1})$ — a clean, deterministic object — yet the value deficit it drives is irreducibly random because the spacings it sums are random. Inconsistency here is not a vague "the bootstrap is off"; it is the precise statement that a deterministic combinatorial law gets composed with a random metric, and no fixed normalization undoes the composition.

The constant is the punchline. The mass $1-e^{-1}$ that fixes the atom is the *same* $1-e^{-1}$ that counts the expected fraction of distinct observations in a resample, that defines Efron's (1983) .632 bootstrap error estimate, and that sets the support-shrinkage rate of iterated subsampling. All three are one fact — $\mathbb{P}(\text{a given point is omitted})=(1-1/n)^n\to e^{-1}$ — viewed from different angles. Here it is the probability the *maximum* is omitted, which is why the bootstrap's chance of even *reaching* $M_n$ is exactly $0.632$, and why the percentile interval, which can never exceed $M_n$, has zero coverage rather than merely poor coverage.

**Limitations.** (i) The exact constants ($\mathrm{Exp}$ spacings, the $(k/n)^n$ identity) are special to a uniform-type endpoint with a positive, finite density at the boundary; for a density vanishing or exploding at the endpoint the index/atom argument (part 1) is unchanged — the rank deficit is *always* Geometric$(1-e^{-1})$, a model-free fact — but the value law (3) rescales by the corresponding extreme-value normalization rather than $1/n$, so $G$ in (4) is endpoint-specific. (ii) Part 3 describes the *average* law; the operative object for a single dataset is the random law of part 2, and our coverage numbers reflect that randomness rather than $G$. (iii) The $m$-out-of-$n$ repair needs a choice of $m$; we exhibit a working band around $\sqrt n$ but do not derive a data-driven optimum (Bickel & Sakov, 2008, give one). (iv) We treat the simplest one-parameter endpoint problem to keep everything exact; the qualitative phenomenon — a counting atom of mass $1-e^{-1}$ at "no improvement" — recurs for any statistic that is a function of a few extreme order statistics.

## References

1. P. J. Bickel and D. A. Freedman, "Some asymptotic theory for the bootstrap," *Annals of Statistics* **9**(6):1196–1217, 1981. DOI:10.1214/aos/1176345637.
2. B. Efron, "Estimating the error rate of a prediction rule: improvement on cross-validation," *Journal of the American Statistical Association* **78**(382):316–331, 1983. DOI:10.1080/01621459.1983.10477973.
3. P. J. Bickel, F. Götze, and W. R. van Zwet, "Resampling fewer than $n$ observations: gains, losses, and remedies for losses," *Statistica Sinica* **7**:1–31, 1997.
4. K. B. Athreya, "Bootstrap of the mean in the infinite variance case," *Annals of Statistics* **15**(2):724–731, 1987. DOI:10.1214/aos/1176350371.
5. D. N. Politis and J. P. Romano, "Large sample confidence regions based on subsamples under minimal assumptions," *Annals of Statistics* **22**(4):2031–2050, 1994. DOI:10.1214/aos/1176325770.
6. P. J. Bickel and A. Sakov, "On the choice of $m$ in the $m$ out of $n$ bootstrap and confidence bounds for extrema," *Statistica Sinica* **18**:967–985, 2008.
7. A. Rényi, "On the theory of order statistics," *Acta Mathematica Academiae Scientiarum Hungaricae* **4**:191–231, 1953. DOI:10.1007/BF02127580.

# Hypergeometric distribution {#hypergeometric_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Hypergeometric distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Sampling without replacement</span>
  <span role="listitem"><strong>C++ type</strong><code>HyperGeometricRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/hypergeometric_pmf.png" loading="lazy" alt="Representative Hypergeometric distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/hypergeometric_cdf.png" loading="lazy" alt="Representative Hypergeometric distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`HyperGeometricRand(N, n, K)` counts successes in \f$n\f$ draws without replacement from a
population of size \f$N\f$ containing \f$K\f$ successes. RandLib requires
\f$N>0\f$, \f$1\leq n\leq N\f$, and \f$1\leq K\leq N\f$. Its support is

\f[
 L=\max(0,n+K-N)\leq k\leq U=\min(n,K),
\f]

and its PMF is

\f[
 P(X=k)=\frac{\binom Kk\binom{N-K}{n-k}}{\binom Nn}.
\f]

## Numerical evaluation and quantiles

The implementation evaluates the log PMF with log-factorials. For \f$F(k)\f$, it compares the
number of support points below and above \f$k\f$, then either sums forward from the true lower
support limit \f$L\f$ or subtracts the shorter upper sum from one. The generating functions are the
finite sums

\f[
 G(z)=\frac1{\binom Nn}
 \sum_{j=L}^{U}\binom Kj\binom{N-K}{n-j}z^j,
 \qquad
 \phi(t)=G(e^{it}).
\f]

The generic finite-support characteristic-function routine evaluates this sum numerically. The

## Moments and shape

principal moments are

\f[
 \operatorname{E}X=\frac{nK}{N},
 \qquad
 \operatorname{Var}X=\frac{nK(N-K)(N-n)}{N^2(N-1)},
\f]

\f[
 \operatorname{Skew}X=
 \frac{(N-2K)(N-2n)\sqrt{N-1}}
 {(N-2)\sqrt{nK(N-K)(N-n)}}.
\f]

Writing \f$A=nK(N-K)(N-n)\f$, the implemented excess kurtosis is

\f[
 \operatorname{ExKurt}X=
 \frac{
  N^2(N-1)\{N(N+1)-6[K(N-K)+n(N-n)]\}
  +6A(5N-6)}
 {A(N-2)(N-3)}.
\f]

The displayed skewness formula applies to nondegenerate populations with \f$N>2\f$, and the
displayed kurtosis formula applies when \f$N>3\f$. For the nondegenerate \f$N=2\f$ and \f$N=3\f$
cases, RandLib evaluates the relevant finite central moment directly; standardized moments are
undefined in deterministic cases. The reported mode is

\f[
 \left\lfloor\frac{(n+1)(K+1)}{N+2}\right\rfloor,
\f]

and the median uses the generic discrete quantile search.

## Sampling

Sampling performs the draws sequentially. After \f$i\f$ objects have been considered and \f$s\f$
successes observed, the next Bernoulli probability is \f$(K-s)/(N-i)\f$; the routine stops early
once every population success has already been selected. Thus the algorithm is exact and uses at
most \f$n\f$ engine draws. There is currently no fitting specialization or dedicated entropy
member.

# Negative-hypergeometric distribution {#negative_hypergeometric_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Negative-hypergeometric distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Sampling without replacement</span>
  <span role="listitem"><strong>C++ type</strong><code>NegativeHyperGeometricRand</code></span>
</div>

## Definition and properties

`NegativeHyperGeometricRand(N, M, m)` counts failures observed before the \f$m\f$th success when

## Sampling

sampling without replacement from \f$N\f$ objects containing \f$M\f$ successes. The constraints
are \f$N>0\f$ and \f$1\leq m\leq M\leq N\f$; the support is
\f$k\in\{0,\ldots,N-M\}\f$. Its PMF is

\f[
 P(X=k)=
 \frac{\binom{k+m-1}{k}\binom{N-m-k}{M-m}}{\binom NM}.
\f]

As in the ordinary hypergeometric implementation, the log PMF uses log-factorials and the CDF
chooses between a forward sum and subtracting the shorter upper sum from one. Its generating
functions are the finite sums

\f[
 G(z)=\frac1{\binom NM}
 \sum_{j=0}^{N-M}
 \binom{j+m-1}{j}\binom{N-m-j}{M-m}z^j,
 \qquad
 \phi(t)=G(e^{it}),
\f]

## Moments and shape

which the generic finite-support routine evaluates numerically. The implemented moments are

\f[
 \operatorname{E}X=\frac{m(N-M)}{M+1},
\f]

\f[
 \operatorname{Var}X=
 \frac{m(N+1)(N-M)(M-m+1)}{(M+1)^2(M+2)}.
\f]

The median and mode use the generic discrete probability searches; skewness, kurtosis, entropy,
and fitting do not have dedicated implementations.

Sampling updates the remaining-success fraction after each draw and stops at the \f$m\f$th
success. If every failure must occur before that stopping success, it returns \f$N-M\f$ early.
The method is exact and requires at most \f$N\f$ sequential Bernoulli draws.

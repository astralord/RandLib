# Chi-squared distribution {#chi_squared_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Chi-squared distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Gamma family</span>
  <span role="listitem"><strong>C++ type</strong><code>ChiSquaredRand</code></span>
</div>

## Definition and properties

`ChiSquaredRand(degree)` takes a positive integer \f$k\f$ and is

\f[
 X\sim\chi_k^2=\operatorname{Gamma}\!\left(\frac{k}{2},\frac12\right).
\f]

Its density, CDF, survival function, characteristic function, quantiles, and sampler are inherited
from the Gamma implementation. In particular,

\f[
 \phi(t)=(1-2it)^{-k/2},\qquad
 \operatorname{E}X=k,\qquad
 \operatorname{Var}X=2k,
\f]

\f[
 \operatorname{Mode}X=\max\{k-2,0\},\qquad
 \operatorname{Skew}X=\sqrt{\frac8k},\qquad
 \operatorname{ExKurt}X=\frac{12}{k}.
\f]

## Fitting

It is also the distribution of \f$\sum_{j=1}^{k}Z_j^2\f$ for independent standard normals.
`fit<ChiSquaredRand>` treats the degree as discrete: it evaluates the log likelihood over positive
integer degrees up to an adaptive, bounded search limit and returns the maximizer. A fixed degree
bypasses the search.

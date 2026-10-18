# Degenerate distribution {#degenerate_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Degenerate distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Bounded algebraic distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>DegenerateRand</code></span>
</div>

## Definition and properties

`DegenerateRand(value)` represents the point mass \f$\Pr[X=a]=1\f$. It is retained in the continuous
class hierarchy for API composition, but it is singular with respect to Lebesgue measure. Its
`f(a)` and `logf(a)` return positive infinity as a sentinel, while both return zero or negative
infinity, respectively, away from \f$a\f$. Its actual distribution functions are

\f[
 F(x)=\mathbf1_{[a,\infty)}(x),\qquad
 S(x)=\mathbf1_{(-\infty,a)}(x),\qquad
 Q(p)=a,
\f]

and

\f[
 \phi(t)=e^{iat},\qquad
 \mathbb{E}[X]=\operatorname{Median}[X]=\operatorname{Mode}[X]=a,
 \qquad \operatorname{Var}[X]=0,
\f]

\f[
 h[X]=0.
\f]

Skewness and excess kurtosis are undefined and return NaN.

## Sampling

## Fitting

Sampling always returns \f$a\f$.
`fit<DegenerateRand>` supports MLE for a nonempty sample only when all observations are exactly
equal; the common observation becomes the fitted point mass.

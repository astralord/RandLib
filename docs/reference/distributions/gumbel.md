# Gumbel distribution {#gumbel_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Gumbel distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Extreme-value and power-tail families</span>
  <span role="listitem"><strong>C++ type</strong><code>GumbelRand</code></span>
</div>

## Definition and properties

`GumbelRand(location, scale)` uses \f$\mu\in\mathbb R\f$ and \f$\beta>0\f$. With
\f$z=(x-\mu)/\beta\f$,

\f[
 f(x)=\frac1\beta e^{-z-e^{-z}},\qquad
 F(x)=e^{-e^{-z}},\qquad
 S(x)=1-e^{-e^{-z}}.
\f]

The survival function uses `expm1`; upper-tail quantiles use `log1p`.

## Sampling

Sampling and the characteristic
function follow from \f$X=\mu-\beta\log E\f$ for a unit exponential \f$E\f$:

\f[
 Q(p)=\mu-\beta\log[-\log p],\qquad
 \phi(t)=e^{i\mu t}\Gamma(1-i\beta t).
\f]

The class currently evaluates the characteristic function through the common expectation quadrature
rather than a dedicated complex-gamma path. Its summaries are

\f[
 \mathbb{E}[X]=\mu+\beta\gamma_E,\qquad
 \operatorname{Var}[X]=\frac{\pi^2\beta^2}{6},
\f]

\f[
 \operatorname{Median}[X]=\mu-\beta\log(\log2),\qquad
 \operatorname{Mode}[X]=\mu,
\f]

\f[
 \operatorname{Skew}[X]=\frac{12\sqrt6\,\zeta(3)}{\pi^3},\qquad
 \operatorname{ExKurt}[X]=\frac{12}{5},\qquad
 h[X]=1+\gamma_E+\log\beta.
\f]

## Fitting

`fit<GumbelRand>` performs constrained MLE, initialized by matching the sample mean and variance to
the two displayed moments.

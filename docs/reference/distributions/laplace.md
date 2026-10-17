# Laplace distribution {#laplace_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Laplace distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Exponential-tail location families</span>
  <span role="listitem"><strong>C++ type</strong><code>LaplaceRand</code></span>
</div>

## Definition and properties

`LaplaceRand(shift, scale)` is the symmetric case \f$\kappa=1\f$:

\f[
 f(x)=\frac1{2\gamma}e^{-|x-m|/\gamma},
 \qquad
 \phi(t)=\frac{e^{imt}}{1+\gamma^2t^2}.
\f]

It uses the asymmetric family's direct piecewise CDF, survival function, logarithmic quantiles, and
difference-of-exponentials sampler; the single-variate fast path equivalently draws an exponential
magnitude and an independent sign bit. Its summaries reduce to

\f[
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=m,
 \qquad
 \operatorname{Var}X=2\gamma^2,
\f]

\f[
 \operatorname{Skew}X=0,\qquad
 \operatorname{ExKurt}X=3,\qquad
 h(X)=1+\log(2\gamma).
\f]

## Fitting

`fit<LaplaceRand>` performs the closed-form MLE: a sample median estimates \f$m\f$, and the mean
absolute deviation from that median estimates \f$\gamma\f$.

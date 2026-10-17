# Uniform distribution {#uniform_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Uniform distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Beta family</span>
  <span role="listitem"><strong>C++ type</strong><code>UniformRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/uniform_pdf.png" loading="lazy" alt="Representative Uniform distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/uniform_cdf.png" loading="lazy" alt="Representative Uniform distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`UniformRand(a, b)` is `BetaRand(1, 1, a, b)`. On \f$[a,b]\f$,

\f[
 f(x)=\frac1d,\qquad F(x)=\frac{x-a}{d},\qquad S(x)=\frac{b-x}{d},
\f]

with the usual constant continuation outside the support, and

\f[
 \phi(t)=
 \begin{cases}
  \dfrac{e^{itb}-e^{ita}}{itd},&t\ne0,\\
  1,&t=0.
 \end{cases}
\f]

## Moments and shape

Its mean, variance, median, skewness, excess kurtosis, and differential entropy are

\f[
 \operatorname{E}X=\operatorname{Median}X=\frac{a+b}{2},\qquad
 \operatorname{Var}X=\frac{d^2}{12},\qquad
 \operatorname{Skew}X=0,
\f]

\f[
 \operatorname{ExKurt}X=-\frac65,\qquad h(X)=\log d.
\f]

Every support point is a mode; `Mode()` returns the midpoint as a representative value.

## Sampling

Sampling is
the affine transform \f$a+dU\f$ of the caller-engine uniform generator. `StandardVariate`,
`StandardVariateClosed`, and `StandardVariateHalfClosed` expose respectively the open, closed, and
left-closed/right-open unit-interval conventions.

## Fitting

`fit<UniformRand>` supports maximum likelihood and UMVU estimation. With order statistics
\f$x_{(1)}\f$ and \f$x_{(n)}\f$, the MLE is

\f[
 \hat a=x_{(1)},\qquad \hat b=x_{(n)},
\f]

while the unbiased estimators are

\f[
 \widetilde a=\frac{nx_{(1)}-x_{(n)}}{n-1},\qquad
 \widetilde b=\frac{nx_{(n)}-x_{(1)}}{n-1}.
\f]

When \f$a\f$ is fixed and the width \f$\theta=b-a\f$ has prior
\f$\operatorname{Pareto}(\alpha,\sigma)\f$, `fit_bayes` returns the conjugate posterior

\f[
 \theta\mid x\sim\operatorname{Pareto}\!\left(
 \alpha+n,\max\{\sigma,x_{(n)}-a\}\right).
\f]

The posterior mean or mode supplies the Bayesian or MAP point estimate, respectively.

# Levy distribution {#levy_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Levy distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Classical heavy-tailed location families</span>
  <span role="listitem"><strong>C++ type</strong><code>LevyRand</code></span>
</div>

## Definition and properties

`LevyRand(location, scale)` uses \f$\mu\in\mathbb R\f$ and \f$\gamma>0\f$, with support closure
\f$[\mu,\infty)\f$. Writing \f$y=x-\mu\f$ for \f$x>\mu\f$,

\f[
 f(x)=\sqrt{\frac{\gamma}{2\pi}}\,y^{-3/2}
 e^{-\gamma/(2y)},
\f]

\f[
 F(x)=\operatorname{erfc}\sqrt{\frac{\gamma}{2y}},\qquad
 S(x)=\operatorname{erf}\sqrt{\frac{\gamma}{2y}}.
\f]

Both tails are therefore direct. Their inverse-error-function forms give

\f[
 Q(p)=\mu+\frac{\gamma}{2[\operatorname{erfc}^{-1}(p)]^2},
\qquad
 Q(1-p)=\mu+\frac{\gamma}{2[\operatorname{erf}^{-1}(p)]^2}.
\f]

## Sampling

The sampler uses the exact normal transform

\f[
 X=\mu+\frac{\gamma}{Z^2},\qquad Z\sim N(0,1),
\f]

and the characteristic function is

\f[
 \phi(t)=\exp\!\left(i\mu t-\sqrt{-2i\gamma t}\right),
\f]

with the principal square root. The mode is \f$\mu+\gamma/3\f$ and the median is the displayed

## Moments and shape

quantile at \f$p=1/2\f$. The mean and variance diverge, while standardized higher moments are
undefined. Its mathematical entropy is

\f[
 h(X)=\frac12+\log(4\gamma\sqrt\pi)+\frac32\gamma_E,
\f]

although the class does not currently expose an entropy member.

With location fixed, the scale MLE is

\f[
 \hat\gamma=\frac{n}{\sum_i(x_i-\mu)^{-1}}.
\f]

## Fitting

`fit<LevyRand>` uses this closed-form scale at its initial location and then performs constrained
two-parameter MLE; the default initial location lies slightly below the sample minimum.

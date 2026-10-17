# Cauchy distribution {#cauchy_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Cauchy distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Classical heavy-tailed location families</span>
  <span role="listitem"><strong>C++ type</strong><code>CauchyRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/cauchy_pdf.png" loading="lazy" alt="Representative Cauchy distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/cauchy_cdf.png" loading="lazy" alt="Representative Cauchy distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`CauchyRand(location, scale)` uses \f$\mu\in\mathbb R\f$ and \f$\gamma>0\f$:

\f[
 f(x)=\frac{\gamma}{\pi[\gamma^2+(x-\mu)^2]},
\qquad
 F(x)=\frac12+\frac1\pi\arctan\frac{x-\mu}{\gamma},
\f]

\f[
 S(x)=\frac12-\frac1\pi\arctan\frac{x-\mu}{\gamma},
\qquad
 Q(p)=\mu-\frac{\gamma}{\tan(\pi p)}.
\f]

The lower and upper tails switch to reciprocal-argument `atan` identities to avoid cancellation.
The characteristic function and entropy are

\f[
 \phi(t)=e^{i\mu t-\gamma|t|},\qquad
 h(X)=\log(4\pi\gamma).
\f]

## Sampling

Sampling draws a uniform point \f$(U,V)\f$ from the unit disk and returns
\f$\mu+\gamma U/V\f$; this ratio is standard Cauchy. The median and mode are \f$\mu\f$. Its mean,
variance, skewness, and kurtosis do not exist; the inherited stable-family API reports NaN for the
mean, skewness, and excess kurtosis and infinity for the variance.

## Fitting

`fit<CauchyRand>` performs constrained two-parameter MLE, initialized at the sample median with
scale equal to half the interquartile range.

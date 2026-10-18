# Exponentially modified Gaussian distribution {#exponentially_modified_gaussian_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Exponentially modified Gaussian distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Normal family</span>
  <span role="listitem"><strong>C++ type</strong><code>ExponentiallyModifiedGaussianRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/emg_pdf.png" loading="lazy" alt="Representative Exponentially modified Gaussian distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/emg_cdf.png" loading="lazy" alt="Representative Exponentially modified Gaussian distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`ExponentiallyModifiedGaussianRand(location, variance, rate)` represents the independent sum

\f[
 Z=X+Y,\qquad X\sim N(\mu,\sigma^2),\qquad
 Y\sim\operatorname{Exponential}(\lambda).
\f]

For all real \f$z\f$,

\f[
 f(z)=\frac\lambda2
 \exp\!\left[\frac\lambda2(2\mu+\lambda\sigma^2-2z)\right]
 \operatorname{erfc}\!\left(
 \frac{\mu+\lambda\sigma^2-z}{\sigma\sqrt2}\right).
\f]

With \f$w=(z-\mu)/\sigma\f$ and \f$v=\lambda\sigma\f$,

\f[
 F(z)=\Phi(w)-e^{-\lambda(z-\mu)+v^2/2}\Phi(w-v),
\f]

\f[
 S(z)=\Phi(-w)+e^{-\lambda(z-\mu)+v^2/2}\Phi(w-v).
\f]

The separately implemented survival expression avoids subtracting the complete CDF from one. Where
\f$w<v\f$, the exponential term is evaluated as \f$\tfrac12e^{-w^2/2}\operatorname{erfcx}((v-w)/\sqrt2)\f$
and, for \f$\mu+\lambda\sigma^2\ge z\f$, the density as
\f$\tfrac\lambda2\operatorname{erfcx}\bigl((\mu+\lambda\sigma^2-z)/(\sqrt2\sigma)\bigr)e^{-w^2/2}\f$, so an
underflowing complementary error function never multiplies an overflowing exponential. The
characteristic function follows directly from independence:

\f[
 \phi(t)=\frac{\exp(i\mu t-\sigma^2t^2/2)}{1-it/\lambda}.
\f]

## Sampling

Sampling adds the library's normal and exponential variates. The mean, variance, skewness, and
excess kurtosis are

\f[
 \mathbb{E}[Z]=\mu+\frac1\lambda,\qquad
 \operatorname{Var}[Z]=\sigma^2+\frac1{\lambda^2},
\f]

\f[
 \operatorname{Skew}[Z]=
 \frac{2\lambda^{-3}}{(\sigma^2+\lambda^{-2})^{3/2}},
 \qquad
 \operatorname{ExKurt}[Z]=
 \frac{6\lambda^{-4}}{(\sigma^2+\lambda^{-2})^2}.
\f]

## Fitting

No elementary median or mode is known. `Median()` uses the generic tail-aware TOMS 748 quantile
solver, and `Mode()` minimizes the negative log density with the shared Brent search, starting from
the mean. There is currently no `fit<ExponentiallyModifiedGaussianRand>` specialization.

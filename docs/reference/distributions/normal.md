# Normal distribution {#normal_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Normal distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Normal family</span>
  <span role="listitem"><strong>C++ type</strong><code>NormalRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/normal_pdf.png" loading="lazy" alt="Representative Normal distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/normal_cdf.png" loading="lazy" alt="Representative Normal distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`NormalRand(location, variance)` uses mean \f$\mu\in\mathbb R\f$ and variance
\f$\sigma^2>0\f$; `FromMeanStd` is the explicit standard-deviation factory. Its density, tails, and
characteristic function are

\f[
 f(x)=\frac{1}{\sigma\sqrt{2\pi}}
 \exp\!\left[-\frac{(x-\mu)^2}{2\sigma^2}\right],
\f]

\f[
 F(x)=\frac12\operatorname{erfc}\!\left(
 \frac{\mu-x}{\sigma\sqrt2}\right),
 \qquad
 S(x)=\frac12\operatorname{erfc}\!\left(
 \frac{x-\mu}{\sigma\sqrt2}\right),
\f]

\f[
 \phi(t)=\exp\!\left(i\mu t-\frac12\sigma^2t^2\right).
\f]

Both tails are evaluated directly with `erfc`. Quantiles use the corresponding direct inverse-tail
forms

\f[
 Q(p)=\mu-\sigma\sqrt2\,\operatorname{erfc}^{-1}(2p),
 \qquad
 Q(1-p)=\mu+\sigma\sqrt2\,\operatorname{erfc}^{-1}(2p).
\f]

RandLib's inverse error function uses piecewise Chebyshev approximations, with separate central,
intermediate-tail, and extreme-tail coefficient sets.

## Sampling

Sampling uses the shared 256-layer normal
ziggurat and then applies \f$\mu+\sigma Z\f$.

The normal law has

\f[
 \mathbb{E}[X]=\operatorname{Median}[X]=\operatorname{Mode}[X]=\mu,
 \qquad
 \operatorname{Var}[X]=\sigma^2,
\f]

\f[
 \operatorname{Skew}[X]=0,\qquad
 \operatorname{ExKurt}[X]=0,\qquad
 h[X]=\frac12\log(2\pi e\sigma^2).
\f]

## Moments and shape

Raw integer moments are evaluated from

\f[
 \mathbb{E}[X^n]=
 \sum_{\substack{0\leq k\leq n\\k\ \mathrm{even}}}
 {n\choose k}\mu^{n-k}\sigma^k(k-1)!!,
\f]

where the \f$k=0\f$ term is one. In RandLib's stable-law convention the same law is
\f$S_1(2,0,\sigma/\sqrt2,\mu)\f$.

## Fitting

`fit<NormalRand>` supports MLE, method of moments, and UMVU estimation. With both parameters free,

\f[
 \hat\mu=\bar x,\qquad
 \hat\sigma^2_{\mathrm{MLE}}=\frac1n\sum_i(x_i-\bar x)^2,\qquad
 \widetilde{\sigma^2}=\frac1{n-1}\sum_i(x_i-\bar x)^2.
\f]

The legacy member API additionally exposes normal/Student-\f$t\f$ confidence intervals for the mean,
chi-squared confidence intervals for the variance, and the gamma-ratio correction that makes the
sample standard deviation unbiased when the mean is known.

The conjugate models exposed by `fit_bayes` are:

- for known \f$\sigma^2\f$ and prior
  \f$\mu\sim N(\mu_0,1/\tau_0)\f$,
  \f$\tau_n=\tau_0+n/\sigma^2\f$ and
  \f$\mu_n=(\tau_0\mu_0+\sum_i x_i/\sigma^2)/\tau_n\f$;
- for known \f$\mu\f$ and prior
  \f$\sigma^2\sim\operatorname{InvGamma}(\alpha_0,\beta_0)\f$,
  \f$\alpha_n=\alpha_0+n/2\f$ and
  \f$\beta_n=\beta_0+\frac12\sum_i(x_i-\mu)^2\f$;
- with both parameters unknown, a normal--inverse-gamma prior updates to

\f[
 \lambda_n=\lambda_0+n,\qquad
 \mu_n=\frac{\lambda_0\mu_0+n\bar x}{\lambda_n},\qquad
 \alpha_n=\alpha_0+\frac n2,
\f]

\f[
 \beta_n=\beta_0+\frac n2\left[
 s_n^2+\frac{\lambda_0}{\lambda_n}(\bar x-\mu_0)^2
 \right],\qquad
 s_n^2=\frac1n\sum_i(x_i-\bar x)^2.
\f]

The joint prior itself is documented under `NormalInverseGammaRand` in the bivariate reference.

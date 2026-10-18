# Gamma distribution {#gamma_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Gamma distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Gamma family</span>
  <span role="listitem"><strong>C++ type</strong><code>GammaRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/gamma_pdf.png" loading="lazy" alt="Representative Gamma distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/gamma_cdf.png" loading="lazy" alt="Representative Gamma distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`GammaRand(shape, rate)` uses the rate parameterization
\f$\alpha>0,\beta>0\f$. The named factories `FromShapeRate` and `FromShapeScale` make the
rate-versus-scale convention explicit; the scale is \f$\theta=1/\beta\f$. For \f$x\geq0\f$,

\f[
 f(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x},\qquad
 F(x)=P(\alpha,\beta x),\qquad
 S(x)=Q(\alpha,\beta x),
\f]

where \f$P\f$ and \f$Q\f$ are the regularized lower and upper incomplete gamma functions. RandLib
has direct implementations of all four of \f$F,S,\log F,\log S\f$; small tails therefore do not
depend on subtracting from one. The characteristic function is

\f[
 \phi(t)=\left(1-\frac{it}{\beta}\right)^{-\alpha}.
\f]

## Moments and shape

The moments and principal summaries are

\f[
 \mathbb{E}[X]=\frac{\alpha}{\beta},\qquad
 \operatorname{Var}[X]=\frac{\alpha}{\beta^2},\qquad
 \operatorname{Mode}[X]=\frac{\max\{\alpha-1,0\}}{\beta},
\f]

\f[
 \operatorname{Skew}[X]=\frac{2}{\sqrt\alpha},\qquad
 \operatorname{ExKurt}[X]=\frac{6}{\alpha}.
\f]

`GeometricMean()` and `GeometricVariance()` return log-moments:

\f[
 \mathbb{E}[\log X]=\psi(\alpha)-\log\beta,
 \qquad \operatorname{Var}[\log X]=\psi_1(\alpha).
\f]

The differential entropy exposed through the exponential-family API is

\f[
 h[X]=\alpha-\log\beta+\log\Gamma(\alpha)+(1-\alpha)\psi(\alpha).
\f]

Except at \f$\alpha=1\f$, where the exponential formula applies, medians and quantiles are found by
inverting the appropriate incomplete-gamma tail. Initial approximations follow Gil, Segura, and
Temme: a power series seeds small lower-tail probabilities, a Lambert-\f$W_{-1}\f$ asymptotic seeds
large values, and an inverse-complementary-error-function plus Lambert-\f$W\f$ approximation handles
large shapes. The result is refined by the shared bracketed TOMS 748 solver. Upper quantiles use
\f$Q\f$ and \f$\log Q\f$ directly.

## Sampling

The sampler dispatches by shape:

- integer shapes one through three are sums of unit exponential variates;
- \f$\alpha=3/2\f$ is an exponential plus half the square of a standard normal;
- very small shapes use Best's RGS algorithm in a configured object and Ahrens--Dieter GS in the
  static `StandardVariate` path;
- \f$1<\alpha<1.2\f$ uses Fishman's rejection algorithm;
- the remaining region uses Marsaglia--Tsang rejection.

The accepted unit-rate variate is divided by \f$\beta\f$. All rejection algorithms consume the
caller-owned engine and enforce the library's defensive iteration limit.

## Fitting

`fit<GammaRand>` performs maximum likelihood for shape and rate on positive samples. If
\f$\bar x\f$ and \f$\overline{\log x}\f$ are the sample means, the profiled rate is
\f$\hat\beta=\hat\alpha/\bar x\f$ and the remaining equation is

\f[
 \log\hat\alpha-\psi(\hat\alpha)
 =\log\bar x-\overline{\log x}.
\f]

RandLib initializes \f$\hat\alpha\f$ with a closed approximation and solves this equation by
safeguarded Newton iteration using \f$1/\alpha-\psi_1(\alpha)\f$. With a known shape, the rate MLE
is \f$\alpha/\bar x\f$ and its UMVU alternative is
\f$(\alpha-1/n)/\bar x\f$.

For fixed shape \f$\alpha\f$, a Gamma prior
\f$\beta\sim\operatorname{Gamma}(\kappa,\gamma)\f$ is conjugate to the rate:

\f[
 \beta\mid x\sim\operatorname{Gamma}\!\left(
 \alpha n+\kappa,\gamma+\sum_i x_i\right).
\f]

`fit_bayes<ExponentialRand>` exposes the \f$\alpha=1\f$ case. The class also exposes the
two-parameter exponential-family representation with sufficient statistic
\f$(\log X,X)\f$ and natural parameter \f$(\alpha-1,-\beta)\f$.

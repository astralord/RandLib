# Exponential distribution {#exponential_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Exponential distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Gamma family</span>
  <span role="listitem"><strong>C++ type</strong><code>ExponentialRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/exponential_pdf.png" loading="lazy" alt="Representative Exponential distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/exponential_cdf.png" loading="lazy" alt="Representative Exponential distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`ExponentialRand(rate)` is the \f$\alpha=1\f$ Gamma law. The factory `FromScale` accepts its mean
\f$\theta=1/\beta\f$. For \f$x\geq0\f$,

\f[
 f(x)=\beta e^{-\beta x},\qquad
 F(x)=1-e^{-\beta x},\qquad
 S(x)=e^{-\beta x},
\f]

\f[
 Q(p)=-\frac{\log(1-p)}{\beta},\qquad
 \phi(t)=\frac{\beta}{\beta-it}.
\f]

## Numerical evaluation and quantiles

RandLib evaluates \f$F\f$ with `expm1`, evaluates \f$S\f$ directly, and uses `log1p` for the lower
quantile so both tails retain precision.

## Sampling

Sampling scales the shared 256-layer unit-exponential
ziggurat by \f$1/\beta\f$.

## Moments and shape

The mean, variance, median, mode, skewness, excess kurtosis, entropy, and raw integer moments are

\f[
 \mathbb{E}[X]=\frac1\beta,\quad
 \operatorname{Var}[X]=\frac1{\beta^2},\quad
 \operatorname{Median}[X]=\frac{\log2}{\beta},\quad
 \operatorname{Mode}[X]=0,
\f]

\f[
 \operatorname{Skew}[X]=2,\qquad
 \operatorname{ExKurt}[X]=6,\qquad
 h[X]=1-\log\beta,\qquad
 \mathbb{E}[X^r]=\frac{r!}{\beta^r}\quad(r=0,1,\ldots).
\f]

## Fitting

`fit<ExponentialRand>` supports the rate MLE \f$1/\bar x\f$ and its UMVU version
\f$(1-1/n)/\bar x\f$. A Gamma prior is conjugate, as described for the Gamma family above.

# Beta-binomial distribution {#beta_binomial_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Beta-binomial distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Binomial family</span>
  <span role="listitem"><strong>C++ type</strong><code>BetaBinomialRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/beta-binomial_pmf.png" loading="lazy" alt="Representative Beta-binomial distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/beta-binomial_cdf.png" loading="lazy" alt="Representative Beta-binomial distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`BetaBinomialRand(number, shape1, shape2)` uses \f$n\in\mathbb N\f$ and
\f$\alpha,\beta>0\f$. It is the beta mixture

\f[
 P\sim\operatorname{Beta}(\alpha,\beta),
 \qquad
 X\mid P\sim\operatorname{Binomial}(n,P),
\f]

and therefore

\f[
 P(X=k)={n\choose k}
 \frac{B(k+\alpha,n-k+\beta)}{B(\alpha,\beta)},
 \qquad k=0,\ldots,n.
\f]

The log PMF is assembled from `lgamma` and log-factorials. The current CDF explicitly sums the
PMF from zero through \f$k\f$; survival and quantiles use the shared discrete paths. The
probability-generating function is

\f[
 G(z)={}_2F_1(-n,\alpha;\alpha+\beta;1-z),
 \qquad
 \phi(t)=G(e^{it}).
\f]

## Moments and shape

The first two moments are

\f[
 \operatorname{E}X=\frac{n\alpha}{\alpha+\beta},
\f]

\f[
 \operatorname{Var}X=
 \frac{n\alpha\beta(\alpha+\beta+n)}
 {(\alpha+\beta)^2(\alpha+\beta+1)}.
\f]

Its skewness is

\f[
 \operatorname{Skew}X=
 \frac{(\alpha+\beta+2n)(\beta-\alpha)}
 {\alpha+\beta+2}
 \sqrt{\frac{\alpha+\beta+1}
 {n\alpha\beta(\alpha+\beta+n)}}.
\f]

For completeness, with \f$s=\alpha+\beta\f$, the implemented excess kurtosis is

\f[
 \operatorname{ExKurt}X=
 \frac{s^2(s+1)}
 {n\alpha\beta(s+2)(s+3)(s+n)}
 \left\{
 s(s-1+6n)+3\left[
 \alpha\beta(n-2)+2n^2
 -\frac{\alpha\beta n(6-n)}s
 -\frac{6\alpha\beta n^2}{s^2}
 \right]\right\}-3.
\f]

The mode
starts from the beta mixing mode scaled by \f$n\f$ and compares the adjacent log probabilities;
the median and characteristic function use the generic finite-support algorithms.

## Sampling

## Fitting

Sampling draws \f$P\f$ from the beta law and then a conditional binomial variate. When
\f$\alpha=\beta=1\f$, the beta-binomial law is exactly uniform on
\f$\{0,\ldots,n\}\f$, so the implementation samples that integer directly. There is currently no
`fit<BetaBinomialRand>` specialization or dedicated entropy member.

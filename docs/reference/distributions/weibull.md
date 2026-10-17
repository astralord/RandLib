# Weibull distribution {#weibull_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Weibull distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Extreme-value and power-tail families</span>
  <span role="listitem"><strong>C++ type</strong><code>WeibullRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/weibull_pdf.png" loading="lazy" alt="Representative Weibull distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/weibull_cdf.png" loading="lazy" alt="Representative Weibull distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`WeibullRand(scale, shape)` uses \f$\lambda,k>0\f$. For \f$x\geq0\f$,

\f[
 f(x)=\frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}
 e^{-(x/\lambda)^k},\qquad
 F(x)=1-e^{-(x/\lambda)^k},\qquad
 S(x)=e^{-(x/\lambda)^k}.
\f]

RandLib handles the finite, zero, or singular density limit at \f$x=0\f$ explicitly. The lower
probability and quantile use `expm1` and `log1p`:

\f[
 Q(p)=\lambda[-\log(1-p)]^{1/k},\qquad
 X=\lambda E^{1/k},\quad E\sim\operatorname{Exponential}(1).
\f]

## Moments and shape

Writing \f$g_r=\Gamma(1+r/k)\f$, the raw moments are
\f$\operatorname{E}X^r=\lambda^r g_r\f$. In particular,

\f[
 \operatorname{E}X=\lambda g_1,\qquad
 \operatorname{Var}X=\lambda^2(g_2-g_1^2),
\f]

\f[
 \operatorname{Median}X=\lambda(\log2)^{1/k},\qquad
 \operatorname{Mode}X=
 \begin{cases}
  0,&k\leq1,\\
  \lambda(1-1/k)^{1/k},&k>1,
 \end{cases}
\f]

and skewness and excess kurtosis are formed from \f$g_1,\ldots,g_4\f$ using the standard central-
moment identities. Its entropy is

\f[
 h(X)=1+\gamma_E(1-1/k)+\log(\lambda/k).
\f]

For \f$k>1\f$, the characteristic function has the convergent power series

\f[
 \phi(t)=\sum_{n=0}^{\infty}
 \frac{\Gamma(1+n/k)}{n!}(i\lambda t)^n.
\f]

RandLib sums its even and odd parts in the log domain when \f$k\geq1\f$ and
\f$\lambda|t|\leq1/2\f$; otherwise it uses expectation quadrature. For \f$k<1\f$, where the
moment series is not convergent, the real integral is split to regularize the density singularity at
zero and the imaginary part is integrated directly.

## Fitting

`fit<WeibullRand>` performs constrained two-parameter MLE, initialized at the exponential model
\f$(\lambda,k)=(\bar x,1)\f$. With \f$k\f$ fixed, the MLE is

\f[
 \hat\lambda=\left(\frac1n\sum_i x_i^k\right)^{1/k}.
\f]

If \f$\lambda^k\sim\operatorname{InverseGamma}(a,b)\f$, the conditional conjugate posterior is

\f[
 \lambda^k\mid x\sim\operatorname{InverseGamma}
 \left(a+n,\ b+\sum_i x_i^k\right),
\f]

which is exposed by `fit_bayes` with the known shape supplied separately.

# Inverse-Gaussian distribution {#inverse_gaussian_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Inverse-Gaussian distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Classical heavy-tailed location families</span>
  <span role="listitem"><strong>C++ type</strong><code>InverseGaussianRand</code></span>
</div>

## Definition and properties

`InverseGaussianRand(mean, shape)` takes \f$\mu,\lambda>0\f$. For \f$x>0\f$,

\f[
 f(x)=\sqrt{\frac{\lambda}{2\pi x^3}}
 \exp\!\left[-\frac{\lambda(x-\mu)^2}{2\mu^2x}\right].
\f]

With

\f[
 a=\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\mu}-1\right),
 \qquad
 b=-\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\mu}+1\right),
\f]

the CDF is

\f[
 F(x)=\Phi(a)+e^{2\lambda/\mu}\Phi(b).
\f]

The naive formula can overflow or lose a tail to subtraction. RandLib combines the common
exponential factor with the scaled complementary error function `erfcx`: it evaluates the lower
tail directly when \f$x<\mu\f$ and the upper tail directly when \f$x\geq\mu\f$.

The characteristic function is

\f[
 \phi(t)=\exp\!\left\{
 \frac{\lambda}{\mu}\left[1-\sqrt{1-\frac{2i\mu^2t}{\lambda}}\right]
 \right\}.
\f]

## Sampling

Sampling uses the Michael--Schucany--Haas transformation. For \f$V=Z^2\f$ with
\f$Z\sim N(0,1)\f$, first form

\f[
 Y=\mu+\frac{\mu^2V}{2\lambda}
 -\frac{\mu}{2\lambda}\sqrt{4\mu\lambda V+\mu^2V^2}.
\f]

## Moments and shape

Return \f$Y\f$ with probability \f$\mu/(\mu+Y)\f$, and otherwise return \f$\mu^2/Y\f$.
The summaries are

\f[
 \operatorname{E}X=\mu,\qquad
 \operatorname{Var}X=\frac{\mu^3}{\lambda},
\f]

\f[
 \operatorname{Mode}X=\mu\left[
 \sqrt{1+\frac{9\mu^2}{4\lambda^2}}-\frac{3\mu}{2\lambda}
 \right],
\f]

\f[
 \operatorname{Skew}X=3\sqrt{\frac{\mu}{\lambda}},\qquad
 \operatorname{ExKurt}X=\frac{15\mu}{\lambda}.
\f]

## Fitting

The median and other quantiles use TOMS 748 over the direct probability functions. The class does
not currently expose a dedicated entropy member. `fit<InverseGaussianRand>` performs constrained
two-parameter MLE, initialized by the sample mean and the moment estimate
\f$\lambda=\bar x^3/s_n^2\f$.

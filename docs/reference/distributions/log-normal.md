# Log-normal distribution {#log_normal_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Log-normal distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Normal family</span>
  <span role="listitem"><strong>C++ type</strong><code>LogNormalRand</code></span>
</div>

## Definition and properties

`LogNormalRand(location, squaredScale)` parameterizes the logarithm:

\f[
 Y=\log X\sim N(\mu,\sigma^2).
\f]

The second constructor argument is the log-variance \f$\sigma^2\f$, while `GetScale()` returns the
log-standard deviation \f$\sigma\f$. The explicit factories `FromLogMeanLogStd` and
`FromLogMeanLogVariance` avoid that distinction at call sites.

For \f$x>0\f$,

\f[
 f(x)=\frac{1}{x\sigma\sqrt{2\pi}}
 \exp\!\left[-\frac{(\log x-\mu)^2}{2\sigma^2}\right],
 \qquad
 F(x)=\Phi\!\left(\frac{\log x-\mu}{\sigma}\right),
\f]

\f[
 S(x)=\Phi\!\left(\frac{\mu-\log x}{\sigma}\right),\qquad
 Q(p)=\exp\!\left[\mu+\sigma\Phi^{-1}(p)\right].
\f]

Probability functions and both quantile tails delegate to the direct normal kernels after taking a
logarithm; sampling exponentiates a normal variate. There is no elementary characteristic function:

\f[
 \phi(t)=\frac{1}{\sigma\sqrt{2\pi}}
 \int_{-\infty}^{\infty}
 \exp\!\left(ite^y-\frac{(y-\mu)^2}{2\sigma^2}\right)dy.
\f]

For \f$\sigma^2<4\f$, RandLib evaluates the real and imaginary expectations with the ordinary
adaptive infinite-interval quadrature. For \f$\sigma^2\geq4\f$, it switches to the half-period
oscillatory integrator with Wynn epsilon acceleration.

Writing \f$v=e^{\sigma^2}\f$, the summaries are

\f[
 \mathbb{E}[X]=e^{\mu+\sigma^2/2},\qquad
 \operatorname{Var}[X]=(v-1)e^{2\mu+\sigma^2},
\f]

\f[
 \operatorname{Median}[X]=e^\mu,\qquad
 \operatorname{Mode}[X]=e^{\mu-\sigma^2},\qquad
 \operatorname{Skew}[X]=(v+2)\sqrt{v-1},
\f]

\f[
 \operatorname{ExKurt}[X]=v^4+2v^3+3v^2-6.
\f]

## Fitting

`fit<LogNormalRand>` is the closed-form normal MLE applied to \f$\log x_i\f$: the fitted location is
the mean log and the fitted squared scale is the population variance of the logs. The same normal,
inverse-gamma, and normal--inverse-gamma conjugate models are available through `fit_bayes`, again
applied in log space.

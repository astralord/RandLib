# Wrapped exponential distribution {#wrapped_exponential_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="other_distribution_reference.html">Circular and singular distributions</a>
  <span aria-hidden="true">/</span>
  <span>Wrapped exponential distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Circular or singular</span>
  <span role="listitem"><strong>Family</strong>Circular distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>WrappedExponentialRand</code></span>
</div>

## Definition and properties

`WrappedExponentialRand(lambda)` uses rate \f$\lambda>0\f$ and represents the exponential law
modulo \f$2\pi\f$ on \f$[0,2\pi]\f$:

\f[
 f(x)=\frac{\lambda e^{-\lambda x}}{1-e^{-2\pi\lambda}},
\f]

\f[
 F(x)=\frac{1-e^{-\lambda x}}{1-e^{-2\pi\lambda}},
 \qquad
 S(x)=\frac{e^{-\lambda x}-e^{-2\pi\lambda}}
 {1-e^{-2\pi\lambda}}.
\f]

`log1mexp` and `expm1` retain accuracy when \f$\lambda\f$ or \f$x\f$ is small, and the survival
function is evaluated directly. For arbitrary real \f$t\f$, the characteristic function is

\f[
 \phi(t)=
 \frac{\lambda\{1-e^{-2\pi(\lambda-it)}\}}
 {(\lambda-it)(1-e^{-2\pi\lambda})}.
\f]

At integer \f$t\f$, the wrap factor cancels and this reduces to
\f$\lambda/(\lambda-it)\f$. The circular summaries are

\f[
 \operatorname{CMean}(X)=\frac\pi2-\arctan\lambda,
 \qquad
 \operatorname{CVar}(X)=1-\frac{\lambda}{\sqrt{1+\lambda^2}}.
\f]

The representative mode is zero and its median is

\f[
 \frac{\log2-\log(1+e^{-2\pi\lambda})}{\lambda}.
\f]

Both lower- and upper-tail quantiles use closed logarithmic forms.

## Sampling

Sampling applies the lower
quantile to one uniform variate, so it is exact and rejection-free. There is currently no fitting
specialization or dedicated entropy member.

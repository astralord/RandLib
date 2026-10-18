# Generalized Pareto distribution {#generalized_pareto_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Generalized Pareto distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Extreme-value and power-tail families</span>
  <span role="listitem"><strong>C++ type</strong><code>GeneralizedParetoRand</code></span>
</div>

## Definition and properties

`GeneralizedParetoRand(shape, scale, location)` uses finite \f$\xi,\mu\f$ and \f$\sigma>0\f$.
For \f$z=(x-\mu)/\sigma\geq0\f$ and \f$T=1+\xi z>0\f$,

\f[
 S(x)=T^{-1/\xi},\qquad
 F(x)=1-T^{-1/\xi},\qquad
 f(x)=\frac1\sigma T^{-1-1/\xi}.
\f]

At \f$\xi=0\f$ these reduce continuously to the shifted exponential law. Negative shape gives the
finite upper endpoint \f$\mu-\sigma/\xi\f$; nonnegative shape gives an infinite right tail.
`log1p`, `expm1`, and direct survival evaluation protect the small-shape and far-tail cases.

## Sampling

Sampling and quantiles use

\f[
 X=
 \begin{cases}
  \mu+\sigma E,&\xi=0,\\
  \mu+\dfrac{\sigma}{\xi}(e^{\xi E}-1),&\xi\ne0,
 \end{cases}
 \qquad E\sim\operatorname{Exponential}(1).
\f]

The characteristic function is evaluated as the defining expectation by the common quadrature.

## Moments and shape

The moments and summaries implemented by the class are

\f[
 \mathbb{E}[X]=\mu+\frac{\sigma}{1-\xi}\quad(\xi<1),
\qquad
 \operatorname{Var}[X]=\frac{\sigma^2}{(1-\xi)^2(1-2\xi)}
 \quad(\xi<1/2),
\f]

\f[
 \operatorname{Median}[X]=
 \begin{cases}
  \mu+\sigma\log2,&\xi=0,\\
  \mu+\dfrac{\sigma}{\xi}(2^\xi-1),&\xi\ne0,
 \end{cases}
\f]

\f[
 \operatorname{Skew}[X]=
 \frac{2(1+\xi)\sqrt{1-2\xi}}{1-3\xi}\quad(\xi<1/3),
\f]

\f[
 \operatorname{ExKurt}[X]=
 \frac{3(1-2\xi)(2\xi^2+\xi+3)}{(1-3\xi)(1-4\xi)}-3
 \quad(\xi<1/4),
 \qquad h[X]=1+\xi+\log\sigma.
\f]

The density is decreasing for \f$\xi>-1\f$, so its mode is \f$\mu\f$. At \f$\xi=-1\f$ it is
constant and the API returns \f$\mu\f$ as a representative mode; for \f$\xi<-1\f$ it increases
toward the finite endpoint and the API returns that endpoint.

## Fitting

`fit<GeneralizedParetoRand>` treats its input as nonnegative exceedances, fixes \f$\mu=0\f$, and
performs constrained two-parameter MLE with \f$\xi>-1/2\f$. A moment estimate initializes the
optimizer. The peaks-over-threshold helper in `randlib/fit/tail.hpp` subtracts the chosen threshold
before invoking this fit.

# Irwin--Hall distribution {#irwin_hall_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Irwin--Hall distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Bounded algebraic distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>IrwinHallRand</code></span>
</div>

## Definition and properties

`IrwinHallRand(number)` is the sum of \f$n\geq1\f$ independent standard uniform variables. On
\f$0\leq x\leq n\f$,

\f[
 f(x)=\frac1{(n-1)!}\sum_{k=0}^{\lfloor x\rfloor}
 (-1)^k{n\choose k}(x-k)^{n-1},
\f]

\f[
 F(x)=\frac1{n!}\sum_{k=0}^{\lfloor x\rfloor}
 (-1)^k{n\choose k}(x-k)^n.
\f]

RandLib has direct branches for \f$n=1\f$ and \f$n=2\f$; the general finite sums use a cached
log-factorial and recursively updated binomial coefficients. Since \f$n-X\f$ has the same law, the
alternating sums are only evaluated for \f$x\le n/2\f$, where they cancel least: the upper tail is
\f$S(x)=F(n-x)\f$ and the density \f$f(x)=f(n-x)\f$. The characteristic function is

\f[
 \phi(t)=\left(\frac{e^{it}-1}{it}\right)^n
 =e^{int/2}\left[\frac{\sin(t/2)}{t/2}\right]^n.
\f]

## Sampling

Sampling explicitly sums \f$n\f$ independent uniforms. Other summaries are

\f[
 \mathbb{E}[X]=\operatorname{Median}[X]=\frac n2,\qquad
 \operatorname{Var}[X]=\frac n{12},\qquad
 \operatorname{Skew}[X]=0,\qquad
 \operatorname{ExKurt}[X]=-\frac6{5n}.
\f]

## Fitting

The API returns \f$n/2\f$ as a representative mode; for \f$n=1\f$ the density is constant and the
mode is not unique. Quantiles use TOMS 748, and there is no `fit<IrwinHallRand>` specialization.

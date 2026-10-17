# Kolmogorov distribution {#kolmogorov_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Kolmogorov distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Spectral and goodness-of-fit laws</span>
  <span role="listitem"><strong>C++ type</strong><code>KolmogorovSmirnovRand</code></span>
</div>

## Definition and properties

`KolmogorovSmirnovRand` is the parameter-free limiting distribution of the two-sided
Kolmogorov--Smirnov statistic. Its CDF has the Jacobi-transformed small-argument series

\f[
 F(x)=\frac{\sqrt{2\pi}}{x}
 \sum_{k=1}^{\infty}
 \exp\!\left[-\frac{(2k-1)^2\pi^2}{8x^2}\right],
 \qquad x>0,
\f]

while its survival function has the rapidly convergent large-argument series

\f[
 S(x)=2\sum_{k=1}^{\infty}(-1)^{k-1}e^{-2k^2x^2}.
\f]

RandLib uses the first representation through \f$x=1\f$ and the second above one, always evaluating
the naturally small tail. Their differentiated density series are

\f[
 f(x)=\frac{\sqrt{2\pi}}{4x^4}
 \sum_{k=1}^{\infty}
 \left[(2k-1)^2\pi^2-4x^2\right]
 e^{-(2k-1)^2\pi^2/(8x^2)}
\f]

for \f$x<1\f$, and

\f[
 f(x)=8x\sum_{k=1}^{\infty}(-1)^{k-1}k^2e^{-2k^2x^2}
\f]

for \f$x\geq1\f$. Every series stops when its next term is negligible relative to the accumulated
sum. `logF` and `logS` use `log1p` only on the complementary side of the switch.

## Moments and shape

The principal summaries are

\f[
 \operatorname{E}X=\sqrt{\frac\pi2}\log2,
 \qquad
 \operatorname{Var}X=\frac{\pi^2}{12}-\frac\pi2(\log2)^2,
\f]

\f[
 \operatorname{Mode}X\approx0.735467812776958,
 \qquad
 \operatorname{Median}X\approx0.82757355518990761.
\f]

The characteristic function and higher standardized moments use the common expectation quadrature.

## Numerical evaluation and quantiles

Quantiles use TOMS 748, initialized by the leading large-argument survival approximation.

## Sampling

Sampling implements Devroye's alternating-series method. It chooses the interval at
\f$x=0.75\f$ using the precomputed probability
\f$F(0.75)\approx0.3728329582237386\f$. The left proposal transforms a rejection sample from a
truncated Gamma \f$(3/2)\f$ law; the right proposal is an exponential transform. Alternating bounds
then accept or reject without evaluating an infinite series to a fixed cutoff. There are no
parameters to fit and no dedicated entropy member.

# Rice distribution {#rice_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Rice distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Noncentral Gaussian norms</span>
  <span role="listitem"><strong>C++ type</strong><code>RiceRand</code></span>
</div>

## Definition and properties

`RiceRand(offset, scale)` uses \f$\nu\geq0\f$ and \f$\sigma>0\f$. It is the magnitude of a
two-dimensional normal vector,

\f[
 X=\sqrt{(\nu+\sigma Z_1)^2+(\sigma Z_2)^2},
 \qquad Z_1,Z_2\stackrel{\mathrm{iid}}{\sim}N(0,1),
\f]

so that

\f[
 \frac{X^2}{\sigma^2}\sim\chi'^2_2\!\left(\frac{\nu^2}{\sigma^2}\right).
\f]

For \f$x>0\f$ its density and direct tails are

\f[
 f(x)=\frac{x}{\sigma^2}
 \exp\!\left[-\frac{x^2+\nu^2}{2\sigma^2}\right]
 I_0\!\left(\frac{x\nu}{\sigma^2}\right),
\f]

\f[
 F(x)=P_1(\nu/\sigma,x/\sigma),\qquad
 S(x)=Q_1(\nu/\sigma,x/\sigma).
\f]

The probability functions use the same saddle-centered Marcum kernel as the noncentral chi-squared
law. When \f$\nu=0\f$ they reduce directly to the Rayleigh `expm1` and exponential formulas. For
the log density, RandLib combines the Bessel exponential with
\f$-(x-\nu)^2/(2\sigma^2)\f$; above Bessel argument 50 it uses the scaled asymptotic series

\f[
 e^{-z}I_0(z)\sim\frac1{\sqrt{2\pi z}}
 \left(1+\frac1{8z}+\frac9{128z^2}
 +\frac{225}{3072z^3}+\frac{11025}{98304z^4}\right),
\f]

which avoids overflow from forming \f$I_0(z)\f$ itself.

## Moments and shape

Let \f$\kappa=\nu^2/(2\sigma^2)\f$ and \f$h=\kappa/2\f$. The implemented first two moments are

\f[
 \mathbb{E}[X]=\sigma\sqrt{\frac\pi2}e^{-h}
 \left[(1+\kappa)I_0(h)+\kappa I_1(h)\right],
\f]

\f[
 \operatorname{Var}[X]=2\sigma^2+\nu^2-(\mathbb{E}[X])^2.
\f]

## Fitting

The mean is evaluated with logarithmic Bessel values and an \f$I_1/I_0\f$ ratio. The sampler is
the defining `hypot` transform. At \f$\nu=0\f$ the mode is \f$\sigma\f$; otherwise the generic
log-density minimizer finds it. Medians and quantiles use TOMS 748, the characteristic function uses
the semi-infinite oscillatory integrator, and higher standardized moments use the common expectation
quadrature. There is currently no `fit<RiceRand>` specialization or dedicated entropy member.

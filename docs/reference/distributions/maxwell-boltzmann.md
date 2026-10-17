# Maxwell--Boltzmann distribution {#maxwell_boltzmann_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Maxwell--Boltzmann distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Nakagami family</span>
  <span role="listitem"><strong>C++ type</strong><code>MaxwellBoltzmannRand</code></span>
</div>

## Definition and properties

`MaxwellBoltzmannRand(scale)` takes \f$\sigma>0\f$ and represents

\f[
 \frac{X}{\sigma}\sim\chi_3,
 \qquad
 X\sim\operatorname{Nakagami}\!\left(\frac32,3\sigma^2\right).
\f]

The factor three in the Nakagami spread is required by the current parameterization and corrects
the missing factor in the archived document. RandLib specializes the density and tails as

\f[
 f(x)=\sqrt{\frac2\pi}\frac{x^2}{\sigma^3}e^{-x^2/(2\sigma^2)},
\f]

\f[
 F(x)=\operatorname{erf}\!\left(\frac{x}{\sigma\sqrt2}\right)
 -\sqrt{\frac2\pi}\frac{x}{\sigma}e^{-x^2/(2\sigma^2)},
\f]

\f[
 S(x)=\operatorname{erfc}\!\left(\frac{x}{\sigma\sqrt2}\right)
 +\sqrt{\frac2\pi}\frac{x}{\sigma}e^{-x^2/(2\sigma^2)}.
\f]

## Moments and shape

The summaries are

\f[
 \operatorname{E}X=2\sigma\sqrt{\frac2\pi},\qquad
 \operatorname{Var}X=\sigma^2\left(3-\frac8\pi\right),\qquad
 \operatorname{Mode}X=\sqrt2\,\sigma,
\f]

\f[
 \operatorname{Skew}X=
 \frac{2\sqrt2(16-5\pi)}{(3\pi-8)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=
 \frac{4\left[\pi(40-3\pi)-96\right]}{(3\pi-8)^2}.
\f]

## Sampling

## Fitting

Sampling uses \f$\sigma\sqrt{2E+Z^2}\f$ with an independent unit exponential \f$E\f$ and standard
normal \f$Z\f$. Quantiles and the characteristic function use the general Nakagami paths. There is
currently no `fit<MaxwellBoltzmannRand>` specialization.

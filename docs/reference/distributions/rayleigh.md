# Rayleigh distribution {#rayleigh_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Rayleigh distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Nakagami family</span>
  <span role="listitem"><strong>C++ type</strong><code>RayleighRand</code></span>
</div>

## Definition and properties

`RayleighRand(scale)` takes \f$\sigma>0\f$ and is

\f[
 \frac{X}{\sigma}\sim\chi_2,
 \qquad
 X\sim\operatorname{Nakagami}(1,2\sigma^2).
\f]

Its direct probability and quantile formulas are

\f[
 f(x)=\frac{x}{\sigma^2}e^{-x^2/(2\sigma^2)},\qquad
 F(x)=1-e^{-x^2/(2\sigma^2)},\qquad
 S(x)=e^{-x^2/(2\sigma^2)},
\f]

\f[
 Q(p)=\sigma\sqrt{-2\log(1-p)},\qquad
 Q(1-p)=\sigma\sqrt{-2\log p}.
\f]

The lower probability uses `expm1` and the lower quantile uses `log1p`.

## Sampling

Sampling uses
\f$\sigma\sqrt{2E}\f$ with the shared unit-exponential ziggurat. The characteristic function is the
\f$\mu=1,\omega=2\sigma^2\f$ case of the Nakagami formula and is evaluated by its quadrature path.

## Moments and shape

The summaries are

\f[
 \operatorname{E}X=\sigma\sqrt{\frac\pi2},\qquad
 \operatorname{Var}X=\frac{4-\pi}{2}\sigma^2,\qquad
 \operatorname{Median}X=\sigma\sqrt{2\log2},\qquad
 \operatorname{Mode}X=\sigma,
\f]

\f[
 \operatorname{Skew}X=
 \frac{2\sqrt\pi(\pi-3)}{(4-\pi)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=
 -\frac{6\pi^2-24\pi+16}{(4-\pi)^2}.
\f]

## Fitting

`fit<RayleighRand>` supports MLE and the unbiased scale estimator. If
\f$R^2=(2n)^{-1}\sum_i x_i^2\f$, then

\f[
 \hat\sigma_{\mathrm{MLE}}=R,\qquad
 \widetilde\sigma=
 \frac{\sqrt n\,\Gamma(n)}{\Gamma(n+1/2)}R.
\f]

The correction follows from
\f$\operatorname{E}R=\sigma\Gamma(n+1/2)/(\sqrt n\,\Gamma(n))\f$. RandLib evaluates the exact
factor in the log domain through \f$n=100\f$ and uses its log-gamma-ratio asymptotic series for
larger samples.

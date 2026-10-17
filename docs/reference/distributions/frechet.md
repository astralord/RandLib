# Frechet distribution {#frechet_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Frechet distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Extreme-value and power-tail families</span>
  <span role="listitem"><strong>C++ type</strong><code>FrechetRand</code></span>
</div>

## Definition and properties

`FrechetRand(shape, scale, location)` uses \f$\alpha,s>0\f$ and \f$m\in\mathbb R\f$, with support
\f$(m,\infty)\f$. For \f$y=(x-m)/s>0\f$,

\f[
 f(x)=\frac{\alpha}{s}y^{-1-\alpha}e^{-y^{-\alpha}},\qquad
 F(x)=e^{-y^{-\alpha}},\qquad
 S(x)=1-e^{-y^{-\alpha}}.
\f]

Survival probabilities again use `expm1`. Direct inversion and sampling give

\f[
 Q(p)=m+s[-\log p]^{-1/\alpha},\qquad
 X=m+sE^{-1/\alpha},\quad E\sim\operatorname{Exponential}(1).
\f]

There is no simpler generally useful characteristic-function expression; RandLib uses its
semi-infinite oscillatory quadrature. Let \f$g_r=\Gamma(1-r/\alpha)\f$. Then

\f[
 \operatorname{E}X=m+sg_1\quad(\alpha>1),\qquad
 \operatorname{Var}X=s^2(g_2-g_1^2)\quad(\alpha>2),
\f]

\f[
 \operatorname{Median}X=m+s(\log2)^{-1/\alpha},\qquad
 \operatorname{Mode}X=m+s\left(\frac{\alpha}{1+\alpha}\right)^{1/\alpha}.
\f]

For \f$\alpha>3\f$ and \f$\alpha>4\f$, respectively,

\f[
 \operatorname{Skew}X=
 \frac{g_3-3g_1g_2+2g_1^3}{(g_2-g_1^2)^{3/2}},
\f]

\f[
 \operatorname{ExKurt}X=
 \frac{g_4-4g_1g_3+6g_1^2g_2-3g_1^4}{(g_2-g_1^2)^2}-3,
 \qquad
 h(X)=1+\gamma_E(1+1/\alpha)+\log(s/\alpha).
\f]

## Fitting

`fit<FrechetRand>` performs constrained three-parameter MLE. Its initializer places the location
slightly below the sample minimum, uses the sample range as scale, and starts with unit shape.

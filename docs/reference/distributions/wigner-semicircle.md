# Wigner semicircle distribution {#wigner_semicircle_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Wigner semicircle distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Bounded algebraic distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>WignerSemicircleRand</code></span>
</div>

## Definition and properties

`WignerSemicircleRand(radius)` takes \f$R>0\f$ and has support \f$[-R,R]\f$:

\f[
 f(x)=\frac{2}{\pi R^2}\sqrt{R^2-x^2},
\f]

\f[
 F(x)=\frac12+\frac1\pi\left[
 \frac{x\sqrt{R^2-x^2}}{R^2}+\arcsin\frac{x}{R}
 \right],\qquad -R\leq x\leq R.
\f]

The implementation does not evaluate that lower tail by adding nearly cancelling terms. For
\f$x\leq0\f$ it sets

\f[
 \delta=1+x/R,\qquad \varphi=2\arcsin\sqrt{\delta/2},\qquad
 F(x)=\frac{2\varphi-\sin(2\varphi)}{2\pi},
\f]

using a local series for \f$z-\sin z\f$ at small \f$z\f$; positive arguments and survival
probabilities follow by symmetry. Quantiles use TOMS 748.

## Sampling

Sampling exploits

\f[
 Y\sim\operatorname{Beta}(3/2,3/2),\qquad X=R(2Y-1).
\f]

The closed-form characteristic function and summaries are

\f[
 \phi(t)=\frac{2J_1(Rt)}{Rt},
\f]

\f[
 \mathbb{E}[X]=\operatorname{Median}[X]=\operatorname{Mode}[X]=0,
 \qquad \operatorname{Var}[X]=\frac{R^2}{4},
\f]

\f[
 \operatorname{Skew}[X]=0,\qquad
 \operatorname{ExKurt}[X]=-1,\qquad
 h[X]=\log(\pi R)-\frac12.
\f]

## Fitting

The class currently obtains \f$\phi\f$ through the common finite-support expectation quadrature
rather than a dedicated Bessel-function path. There is no `fit<WignerSemicircleRand>`
specialization.

# Arcsine subfamily {#arcsine_subfamily_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Arcsine subfamily</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Beta family</span>
  <span role="listitem"><strong>C++ type</strong><code>ArcsineRand</code></span>
</div>

## Definition and properties

`ArcsineRand(shape, a, b)` accepts \f$0<q<1\f$ and is the Beta subfamily

\f[
 X\sim\operatorname{Beta}(1-q,q,a,b).
\f]

Thus, for \f$0<y<1\f$,

\f[
 f(x)=\frac{\sin(\pi q)}{\pi d}\,y^{-q}(1-y)^{q-1},\qquad
 F(x)=I_y(1-q,q),
\f]

\f[
 \phi(t)=e^{ita}\,{}_1F_1(1-q;1;itd),\qquad
 \operatorname{E}X=a+d(1-q),\qquad
 \operatorname{Var}X=\frac{d^2q(1-q)}{2}.
\f]

Only the default \f$q=1/2\f$ is the classical symmetric arcsine distribution. In that case,

\f[
 F(x)=\frac{2}{\pi}\arcsin\sqrt y,\qquad
 S(x)=\frac{2}{\pi}\arccos\sqrt y,\qquad
 Q(p)=a+d\sin^2\!\left(\frac{\pi p}{2}\right),
\f]

and RandLib uses these formulas directly, including the inverse-CDF sampler. Other shapes inherit
the Beta probability, quadrature, quantile, and sampling paths.

## Fitting

`fit<ArcsineRand>` estimates \f$q\f$ by maximum likelihood on the standard support. If
\f$L_0=\overline{\log Y}\f$ and \f$L_1=\overline{\log(1-Y)}\f$, the score equation is

\f[
 L_1-L_0+\pi\cot(\pi q)=0,
\f]

which the implementation solves in closed form with an `atan` expression and the appropriate branch
in \f$(0,1)\f$.

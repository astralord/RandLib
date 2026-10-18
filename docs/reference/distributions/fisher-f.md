# Fisher F distribution {#fisher_f_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Fisher F distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Transformed beta and gamma laws</span>
  <span role="listitem"><strong>C++ type</strong><code>FisherFRand</code></span>
</div>

## Definition and properties

`FisherFRand(degree1, degree2)` requires positive integer degrees \f$d_1,d_2\f$. It is a scaled
beta-prime law,

\f[
 \frac{d_1}{d_2}X\sim
 \operatorname{BetaPrime}\!\left(\frac{d_1}{2},\frac{d_2}{2}\right),
\f]

and therefore

\f[
 f(x)=\frac{(d_1/d_2)^{d_1/2}x^{d_1/2-1}}
 {B(d_1/2,d_2/2)\left(1+d_1x/d_2\right)^{(d_1+d_2)/2}},
\f]

\f[
 F(x)=I_{d_1x/(d_1x+d_2)}\!\left(\frac{d_1}{2},\frac{d_2}{2}\right),
 \qquad
 S(x)=I_{d_2/(d_1x+d_2)}\!\left(\frac{d_2}{2},\frac{d_1}{2}\right).
\f]

Its characteristic function is the beta-prime characteristic function evaluated at
\f$td_2/d_1\f$ and is computed by the same oscillatory integrator. Quantiles and samples are likewise
scaled beta-prime results.

The main summaries are

\f[
 \mathbb{E}[X]=\frac{d_2}{d_2-2}\quad(d_2>2),
\f]

\f[
 \operatorname{Var}[X]=
 \frac{2d_2^2(d_1+d_2-2)}{d_1(d_2-2)^2(d_2-4)}
 \quad(d_2>4),
\f]

\f[
 \operatorname{Mode}[X]=
 \begin{cases}
 \dfrac{d_2(d_1-2)}{d_1(d_2+2)},&d_1>2,\\
 0,&d_1\leq2,
 \end{cases}
\f]

\f[
 \operatorname{Skew}[X]=
 \frac{2d_1+d_2-2}{d_2-6}
 \sqrt{\frac{8(d_2-4)}{d_1(d_1+d_2-2)}}
 \quad(d_2>6),
\f]

\f[
 \operatorname{ExKurt}[X]=
 \frac{12}{(d_2-6)(d_2-8)}
 \left[
 \frac{(d_2-2)^2(d_2-4)}{d_1(d_1+d_2-2)}+5d_2-22
 \right]
 \quad(d_2>8).
\f]

Undefined higher moments are reported as infinity.

## Fitting

`fit<FisherFRand>` treats both degrees as
discrete and searches the likelihood grid from one through 80 for each unfixed degree. This finite
implementation search bound is part of the fitting API behavior, not a property of the law.

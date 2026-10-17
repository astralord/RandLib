# Holtsmark and Landau distributions {#holtsmark_landau_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Holtsmark and Landau distributions</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Stable families</span>
  <span role="listitem"><strong>C++ type</strong><code>HoltsmarkRand</code>, <code>LandauRand</code></span>
</div>

## Definition and properties

The named stable specializations preserve their conventional S1 definitions:

\f[
 \operatorname{Holtsmark}(\gamma,\mu)
 =S_1\!\left(\frac32,0,\gamma,\mu\right),
 \qquad
 \phi(t)=\exp\!\left(i\mu t-|\gamma t|^{3/2}\right),
\f]

\f[
 \operatorname{Landau}(\gamma,\mu)
 =S_1(1,1,\gamma,\mu),
\f]

\f[
 \phi(t)=
 \exp\!\left\{
  i\mu t-\gamma|t|
  \left[1+i\frac{2}{\pi}\operatorname{sgn}(t)\log|t|\right]
 \right\}.
\f]

Holtsmark is symmetric, so its mean, median, and mode are \f$\mu\f$; its variance is infinite.
Landau has full real support, a divergent positive mean and infinite variance; its median and mode
are obtained numerically. Both classes reuse every probability, quantile, and sampling path of
`StableDistribution`. They do not have separate fitting specializations.

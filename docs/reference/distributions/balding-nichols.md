# Balding--Nichols distribution {#balding_nichols_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Balding--Nichols distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Beta family</span>
  <span role="listitem"><strong>C++ type</strong><code>BaldingNicholsRand</code></span>
</div>

## Definition and properties

`BaldingNicholsRand(fixationIndex, frequency)` accepts \f$F,p\in(0,1)\f$ and is the standard-support
Beta law

\f[
 X\sim\operatorname{Beta}\!\left(p\frac{1-F}{F},
 (1-p)\frac{1-F}{F}\right).
\f]

## Fitting

Consequently \f$\operatorname{E}X=p\f$ and \f$\operatorname{Var}X=Fp(1-p)\f$; its density, CDF,
characteristic function, higher moments, quantiles, and sampler are the corresponding Beta formulas
and algorithms above. The constructor order is \f$(F,p)\f$, despite the conventional notation often
listing frequency first. There is currently no `fit<BaldingNicholsRand>` specialization.

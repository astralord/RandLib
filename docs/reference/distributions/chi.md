# Chi distribution {#chi_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Chi distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Nakagami family</span>
  <span role="listitem"><strong>C++ type</strong><code>ChiRand</code></span>
</div>

## Definition and properties

`ChiRand(degree)` takes a positive integer \f$k\f$ and is

\f[
 X\sim\chi_k=\operatorname{Nakagami}\!\left(\frac{k}{2},k\right),
 \qquad X^2\sim\chi_k^2.
\f]

Thus

\f[
 f(x)=\frac{2^{1-k/2}}{\Gamma(k/2)}x^{k-1}e^{-x^2/2},
 \qquad
 F(x)=P\!\left(\frac{k}{2},\frac{x^2}{2}\right),
\f]

\f[
 \operatorname{E}X=\sqrt2\frac{\Gamma((k+1)/2)}{\Gamma(k/2)},\qquad
 \operatorname{Var}X=k-(\operatorname{E}X)^2,\qquad
 \operatorname{Mode}X=\sqrt{k-1}.
\f]

## Fitting

Its higher moments, characteristic function, direct upper tail, quantiles, and sampler are the
corresponding Nakagami operations above. There is currently no `fit<ChiRand>` specialization.

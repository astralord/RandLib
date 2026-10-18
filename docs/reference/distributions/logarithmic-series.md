# Logarithmic-series distribution {#logarithmic_series_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Logarithmic-series distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Geometric and logarithmic-series families</span>
  <span role="listitem"><strong>C++ type</strong><code>LogarithmicRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/logarithmic_pmf.png" loading="lazy" alt="Representative Logarithmic-series distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/logarithmic_cdf.png" loading="lazy" alt="Representative Logarithmic-series distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

For \f$0<p<1\f$ and \f$k\in\{1,2,\ldots\}\f$,

\f[
 \Pr[X=k]=-\frac{p^k}{k\log(1-p)}.
\f]

Writing

\f[
 B_p(a,0)=\int_0^p\frac{t^{a-1}}{1-t}\,dt
         =p^a\sum_{j=0}^{\infty}\frac{p^j}{a+j},
\f]

the implemented tails are

\f[
 F(k)=1+\frac{B_p(k+1,0)}{\log(1-p)},
 \qquad
 S(k)=-\frac{B_p(k+1,0)}{\log(1-p)}.
\f]

The series routine writes out its first three terms and continues until the next term is negligible;
the survival probability is evaluated directly from that positive tail.

The probability-generating and characteristic functions are

\f[
 G(z)=\frac{\log(1-pz)}{\log(1-p)},
 \qquad
 \phi(t)=\frac{\log(1-pe^{it})}{\log(1-p)}.
\f]

## Moments and shape

With \f$q=1-p\f$, the implemented moments are

\f[
 \mathbb{E}[X]=-\frac{p}{q\log q},
 \qquad
 \operatorname{Var}[X]=-\frac{p(p+\log q)}{q^2(\log q)^2}.
\f]

The mode is one and the median is obtained by the generic discrete quantile search.

## Sampling

Sampling uses
Kemp's second accelerated generator: two uniforms select the atom at one, the atom at two, or an
inverse logarithmic tail in constant expected time. There is currently no fitting specialization or
dedicated entropy member.

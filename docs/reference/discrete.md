# Discrete distributions {#discrete_distribution_reference}

The discrete reference covers `BernoulliRand`, `BetaBinomialRand`, `BinomialRand`, `CategoricalRand`,
`GeometricRand`, `HyperGeometricRand`, `LogarithmicRand`, `NegativeBinomialRand`,
`NegativeHyperGeometricRand`, `PoissonRand`, `RademacherRand`, `SkellamRand`, `UniformDiscreteRand`,
`YuleRand`, `ZetaRand` and `ZipfRand`.

For every class, `F(k)` means \f$\Pr[X\leq k]\f$ and `S(k)` means \f$\Pr[X>k]\f$. Classes with a
dedicated `S` evaluate that upper tail directly; otherwise the shared implementation returns
\f$1-F(k)\f$. Discrete quantiles are the smallest integers whose CDF reaches the requested
probability. The generic characteristic-function path sums outward from the median and uses an
Euler-transformed remainder for unbounded slowly decaying tails.

## Binomial family

<div class="reference-callout">
  <strong>16 focused entries</strong>
  <span>Choose a family below; each distribution opens on its own page with definitions, formulas, numerical methods, and fitting notes.</span>
</div>

## Binomial family {#discrete_binomial_family}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="bernoulli_distribution.html"><strong>Bernoulli distribution</strong></a>
    <code>BernoulliRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="binomial_distribution.html"><strong>Binomial distribution</strong></a>
    <code>BinomialRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="beta_binomial_distribution.html"><strong>Beta-binomial distribution</strong></a>
    <code>BetaBinomialRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Finite categorical laws {#discrete_finite_categorical_laws}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="categorical_distribution.html"><strong>Categorical distribution</strong></a>
    <code>CategoricalRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="discrete_uniform_distribution.html"><strong>Discrete-uniform distribution</strong></a>
    <code>UniformDiscreteRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="rademacher_distribution.html"><strong>Rademacher distribution</strong></a>
    <code>RademacherRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Poisson family {#discrete_poisson_family}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="poisson_distribution.html"><strong>Poisson distribution</strong></a>
    <code>PoissonRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="skellam_distribution.html"><strong>Skellam distribution</strong></a>
    <code>SkellamRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Geometric and logarithmic-series families {#discrete_geometric_logarithmic_series_families}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="geometric_distribution.html"><strong>Geometric distribution</strong></a>
    <code>GeometricRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="negative_binomial_distribution.html"><strong>Negative-binomial distribution</strong></a>
    <code>NegativeBinomialRand, PascalRand, PolyaRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="logarithmic_series_distribution.html"><strong>Logarithmic-series distribution</strong></a>
    <code>LogarithmicRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Sampling without replacement {#discrete_sampling_without_replacement}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="hypergeometric_distribution.html"><strong>Hypergeometric distribution</strong></a>
    <code>HyperGeometricRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="negative_hypergeometric_distribution.html"><strong>Negative-hypergeometric distribution</strong></a>
    <code>NegativeHyperGeometricRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Discrete power-law families {#discrete_discrete_power_law_families}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="zeta_distribution.html"><strong>Zeta distribution</strong></a>
    <code>ZetaRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="truncated_zipf_distribution.html"><strong>Truncated Zipf distribution</strong></a>
    <code>ZipfRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Discrete</span>
    <a href="yule_simon_distribution.html"><strong>Yule--Simon distribution</strong></a>
    <code>YuleRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Related guides

- @ref mathematical_reference "Mathematical-reference conventions"
- @ref numerical_methods_reference "Shared numerical methods"
- @ref md_docs_2fitting "Fitting API"

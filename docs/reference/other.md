# Circular and singular distributions {#other_distribution_reference}

This page covers the circular `VonMisesRand` and `WrappedExponentialRand` laws and the singular
`CantorRand` law. Their support and probability operations do not fit the ordinary real-line density
template, so each entry explicitly states its measure and CDF conventions.

## Circular distributions

For circular laws, ordinary linear moments depend on the chosen cut of the circle and RandLib
returns NaN from `Mean()`, `Variance()`, `Skewness()` and `ExcessKurtosis()`. `CircularMean()` is an
angle and `CircularVariance()` is \f$1-|\mathbb{E}[e^{iX}]|\f$. Densities and CDFs use the
single representative interval reported by `MinValue()` and `MaxValue()`.

<div class="reference-callout">
  <strong>3 focused entries</strong>
  <span>Choose a family below; each distribution opens on its own page with definitions, formulas, numerical methods, and fitting notes.</span>
</div>

## Circular distributions {#circular-singular_circular}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Circular or singular</span>
    <a href="von_mises_distribution.html"><strong>Von Mises distribution</strong></a>
    <code>VonMisesRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Circular or singular</span>
    <a href="wrapped_exponential_distribution.html"><strong>Wrapped exponential distribution</strong></a>
    <code>WrappedExponentialRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Singular distributions {#circular-singular_singular}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Circular or singular</span>
    <a href="cantor_distribution.html"><strong>Cantor distribution</strong></a>
    <code>CantorRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Related guides

- @ref mathematical_reference "Mathematical-reference conventions"
- @ref numerical_methods_reference "Shared numerical methods"
- @ref md_docs_2fitting "Fitting API"

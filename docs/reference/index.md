# Mathematical reference {#mathematical_reference}

This reference describes the probability laws implemented by RandLib and the algorithms used to
evaluate and sample them. Its statements are checked against the public headers and tests.

RandLib currently exposes 72 concrete distribution classes: 50 continuous, 16 discrete, three
bivariate, two circular and one singular. Related classes that share one mathematical family are
documented together. Internal base classes are discussed only when their parameterization or
algorithm affects a public law.

<div class="reference-hero"><strong>Find a distribution</strong><span>Browse by kind and then by mathematical family. Every law has a focused page with its formulas and implementation notes.</span></div>

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">50 C++ types · 48 entries</span>
    <a href="continuous_distribution_reference.html"><strong>Continuous distributions</strong></a>
    <span>Densities, tails, quantiles, moments, samplers, and fitting.</span>
    <span class="distribution-card__meta">12 mathematical families</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">16 C++ types · 16 entries</span>
    <a href="discrete_distribution_reference.html"><strong>Discrete distributions</strong></a>
    <span>Mass functions, generating functions, finite and infinite supports.</span>
    <span class="distribution-card__meta">6 mathematical families</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">3 C++ types</span>
    <a href="bivariate_distribution_reference.html"><strong>Bivariate distributions</strong></a>
    <span>Joint densities and masses, lower-orthant probabilities, and covariance.</span>
    <span class="distribution-card__meta">3 focused entries</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">3 C++ types</span>
    <a href="other_distribution_reference.html"><strong>Circular and singular laws</strong></a>
    <span>Periodic densities and distributions without a Lebesgue density.</span>
    <span class="distribution-card__meta">3 focused entries</span>
  </div>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/normal_pdf.png" loading="lazy" alt="Normal density curves"><span class="distribution-plot__caption">Continuous density: normal family.</span></div><div class="distribution-plot"><img src="images/poisson_pmf.png" loading="lazy" alt="Poisson probability mass functions"><span class="distribution-plot__caption">Discrete mass: Poisson family.</span></div><div class="distribution-plot"><img src="images/stable_pdf.png" loading="lazy" alt="Stable density curves"><span class="distribution-plot__caption">Heavy tails and skewness: stable family.</span></div></div>

## How to read an entry

Each distribution entry records, when applicable:

- the RandLib class and constructor parameterization;
- support, probability density or mass, CDF and direct survival function;
- characteristic function or probability-generating function;
- mean, variance and principal higher moments, including existence conditions;
- relations to other implemented laws;
- the sampling algorithm and caller-owned random-engine behavior;
- the numerical methods behind probability functions and quantiles;
- fitting methods and important limitations.

The formulas use TeX notation and are rendered by MathJax in the generated Doxygen site. Source
Markdown retains the TeX directly, so the mathematical content does not depend on generated image
files.

## Notation and API correspondence

For a continuous random variable \f$X\f$,

\f[
 f(x)=\frac{d}{dx}\Pr[X\le x],\qquad
 F(x)=\Pr[X\le x],\qquad
 S(x)=\Pr[X>x],\qquad
 \phi(t)=\mathbb{E}[e^{itX}].
\f]

For a discrete random variable, \f$p(k)=\Pr[X=k]\f$ replaces the density. RandLib uses `f` for a
continuous density, `P` for a discrete mass, `F` for the CDF, `S` for the survival function, `Variate`
or `operator()` for one draw, and `CF` for the characteristic function. `Quantile(p)` is the lower
quantile

\f[
 Q(p)=\inf\{x:F(x)\ge p\}.
\f]

The API deliberately evaluates upper tails directly when subtraction \f$1-F(x)\f$ would lose
precision. Formulas in this reference may state \f$S=1-F\f$ as a mathematical identity even when
the implementation uses a separate series, special function or asymptotic expansion.

Unless an entry says otherwise, `Variance()` returns the second central moment, `Skewness()` the
standardized third central moment and `ExcessKurtosis()` the standardized fourth central moment
minus three. Infinity denotes a divergent moment; NaN denotes a standardized moment that is not
defined.

## Further guides

- @ref numerical_methods_reference "Shared numerical methods"
- @ref md_docs_2stable-distributions "Stable-law parameterizations and accuracy"
- @ref md_docs_2fitting "Fitting, censoring, truncation, and diagnostics"

## Source and evidence policy

The implementation is authoritative for RandLib parameter names and conventions; a standard
mathematical source is authoritative for the law itself. Numerical claims follow the
[accuracy and reference-data policy](../accuracy-policy.md) and are backed by exact identities,
committed high-precision data or focused regression tests.

The focused pages currently publish 54 representative density, mass, and CDF plots migrated from
the archived PDF for the laws whose constructor conventions still match. Plot legends identify the
illustrated parameter sets; the current formulas and constructor descriptions remain authoritative.

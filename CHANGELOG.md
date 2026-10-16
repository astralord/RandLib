# Changelog

All notable changes to RandLib are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [2.1.0] - Unreleased

### Stable distributions

- Added Nolan's S0 parameterization alongside the source-compatible S1 default. Constructors,
  location accessors, setters, names, support endpoints, probability functions, characteristic
  functions and sampling consistently use the selected convention, including the special
  logarithmic location conversion at alpha = 1.
- Implemented stable densities, distribution functions and quantiles throughout the former skewed
  near-alpha-one gap. A cancellation-free S0 characteristic function and alpha-one sampling limit
  keep adjacent representable exponents continuous; 40-digit mpmath values cover both sides.
- Added McCulloch quantile fitting as a fast robust estimator and initializer for constrained stable
  maximum likelihood. Fits are reported in the continuous S0 parameterization.
- Added alpha-one density and probability tail expansions through eighth order. Geometric-stable
  mixture evaluation now covers its full historical tail grid in seconds instead of minutes, with
  focused 40-digit mpmath references on both sides.
- Added the log-stable distribution as the exponential transform of a stable law, with transformed
  density, direct tails, sampling, quantiles, median, correct boundary/interior modes, S0/S1
  parameterization and stable-family fitting. Exact Laplace-exponent moments cover totally
  left-skewed laws and the Gaussian case; other non-Gaussian skewnesses report divergent moments.
  Its characteristic function deliberately throws `not_implemented` away from zero.

### Distributions

- Added the Tier 1 noncentral t distribution with centered Benton–Krishnamoorthy probability
  series, direct survival probabilities, recurrence and integral densities, transformed sampling,
  quantiles, moments, a conditional-normal characteristic function and constrained maximum-
  likelihood fitting. Saddle-split integration and a small-denominator expansion retain accuracy
  for large noncentrality and representable extreme tails; a 243-value 40-digit mpmath table and
  focused exact-series regressions provide evidence.
- Made its centered series and moments cancellation-resistant on platforms where `long double` is
  binary64, using Loader's Poisson saddle-point form, stable log-beta ratios, compensated summation
  and grouped central-moment identities. Exact probability levels at zero now invert to zero.
- Evaluated its moderate-argument density directly from the hypergeometric representation, avoiding
  cancellation between two probabilities on binary64-only platforms.

### Numerical kernel

- Added type-generic tanh-sinh quadrature for finite intervals with integrable endpoint
  singularities. Endpoint-aware callables receive the accurately computed distance to the nearer
  boundary, avoiding cancellation when the abscissa itself rounds to that boundary.
- Added QAWF-style semi-infinite Fourier quadrature with half-period subdivision and Wynn epsilon
  acceleration. Heavy-tailed Pareto, Frechet, inverse-gamma, Fisher F and beta-prime characteristic
  functions, plus the large-variance log-normal case, now meet the standard `1e-8` absolute target;
  all twelve former tolerance overrides were removed.
- Added direct Marcum-Q evaluation through a saddle-centered Poisson mixture. Noncentral chi-squared
  survival probabilities now retain relative accuracy through representable tails instead of
  subtracting the CDF from one; all three former tolerance overrides were removed.
- Added the Tier 2 Rice distribution with direct density, CDF and survival functions, transformed-
  normal sampling, closed-form mean and variance, and 40-digit mpmath reference data.

### Build and tooling

- Avoided a GCC 12 false positive in BFGS identity-matrix construction and replaced `std::jthread`
  in the parallel-stream example for Apple libc++ versions that do not provide it.
- Reduced the stable-law tangent argument before evaluation and made quarter-turn values exact,
  keeping totally skewed short-tail quantiles portable when `long double` is binary64.

## [2.0.0] - Unreleased

### Documentation

- Rewrote the README around distribution fitting and heavy-tail work, with a quick start and an
  explicit comparison with the standard library, Boost.Math and StatsLib.
- Added guides for fitting, parallel stream ownership, stable-law parameterization and limitations,
  the numerical evidence policy, packaging, and every v1-to-v2 breaking change.
- Added a strict Doxygen build with the pinned Doxygen Awesome theme and a GitHub Pages deployment
  workflow. Redundant conversion-era `@fn` directives and incomplete public parameter comments were
  cleaned up so documentation warnings fail the build.
- Retired the generated legacy PDF. Its LaTeX source and figures remain under `docs/legacy/` as a
  clearly marked 1.x mathematical reference.

### Examples and packaging

- Added warning-clean, CI-compiled examples for standard-engine sampling, parallel streams, fitted
  model selection, and conjugate Bayesian updates. The stable-fitting example follows the Phase 5
  stable estimator.
- Added a locally validated vcpkg overlay port and Conan 2 header-only recipe with a compiled
  `test_package` consumer.
- Added deterministic single-header generation and a CI smoke test. Tagged v2 releases extract their
  matching changelog section and attach `randlib_single.hpp` with its SHA-256 checksum.
- Removed a Bernoulli/binomial include cycle exposed by standalone public-header consumption without
  changing the exact sampling algorithm or engine advancement.

## [2.0.0-beta.1] - Unreleased

### Fitting

- Began the unified non-mutating fitting API with `FitResult`, `FitOptions`, likelihood-based
  AIC/BIC, observed-information covariance and standard errors, and fixed-parameter MLE. Existing
  estimators are exposed for normal, log-normal, exponential, gamma, beta, arcsine, beta-prime,
  Pareto, Laplace, Lévy, Rayleigh, uniform, degenerate, Poisson, binomial, negative-binomial and
  discrete-uniform distributions. The generic constrained optimizer supplies full Weibull and
  logistic MLEs, plus MLEs for Student t, Cauchy, Gumbel, Fréchet, inverse-gamma,
  inverse-Gaussian and noncentral chi-squared laws. Bernoulli and geometric fits are closed form;
  chi-squared, Erlang and Fisher F use bounded integer likelihood searches.
- Added named method-of-moments fitting, conjugate Bayesian posterior updates for the existing
  normal, log-normal, exponential, Poisson, binomial, negative-binomial, Weibull, Pareto and uniform
  models, and AIC/BIC model selection across candidate distributions.
- Added continuous goodness-of-fit results for Kolmogorov--Smirnov, Anderson--Darling and
  Cramer--von Mises statistics with p-values, plus a discrete chi-squared test that combines sparse
  bins and accounts for fitted parameters.
- Added seeded parameter-recovery grids over multiple sample sizes and committed fixed datasets with
  independently generated SciPy MLE results for nine continuous distributions. Reference-tool
  versions are pinned in `tools/requirements.txt`.
- Added generalized Pareto and generalized extreme-value distributions, with high-precision mpmath
  reference data, generic MLEs, and stable direct survival-function implementations.
- Added Hill and Pickands tail-index estimators, plot-ready Hill diagnostics, and
  peaks-over-threshold generalized-Pareto fitting.
- Added maximum-likelihood fitting for exact, left-, right- and interval-censored observations and
  for one- or two-sided truncated samples. These adapters share the constrained optimizer,
  information criteria and observed-information covariance used by ordinary MLE.
- Censored and truncated likelihood fits use a simplex search from the original initializer,
  avoiding gradient-based false convergence on flat transformed-likelihood regions.

### API compatibility

- Marked the legacy mutating member `Fit*` methods deprecated in favor of `fit`, `fit_moments` and
  `fit_bayes`. They remain source-compatible until v3; projects temporarily unable to migrate can
  define `RANDLIB_DISABLE_DEPRECATED_FIT_WARNINGS` before including RandLib.
- Removed the unused private `POISSON` binomial generator selector. The exact Bernoulli-sum,
  waiting-time and rejection samplers remain; no approximation is introduced.

## [2.0.0-alpha.1] - Unreleased

### Architecture

- Converted RandLib to a header-only C++20 library under `include/randlib/` and moved the public API
  into namespace `randlib`.
- Made every sampler accept any standard uniform random bit generator. Added PCG64, Xoshiro256++,
  SplitMix64 and deterministic stream construction while retaining the corrected legacy engines.
- Removed hidden mutable sampling state. Bulk sampling uses `std::span`, support categories are a
  scoped enum, parameter factories are named, and distribution classes expose const-correct,
  `[[nodiscard]]` APIs.

### Build and tooling

- Added a pinned Clang 18 formatting policy and pre-commit hook, and applied it to the public headers
  and test sources without reordering dependency-sensitive includes.
- Clang-tidy now analyzes the explicit-instantiation translation unit and public headers as a blocking
  CI job, with every diagnostic promoted to an error.
- Added opt-in Google Benchmark sampling comparisons for normal, exponential, gamma and Poisson
  distributions against the standard library, plus stable-law sampling.
- Engine-adapter interoperability tests now verify identical engine advancement and allow up to
  eight ULPs of platform-dependent variation between direct and type-erased sampling paths.

### API cleanup

- Kept the type-erased `EngineRef` sampling adapters private and removed the access-specifier churn
  left by the automated header conversion.
- Moved implementation-only numerical helpers and tolerances into `randlib::math::detail`.
  The public `Pair` and `Triplet` aliases now live with the other core declarations, while the
  unconstrained arithmetic overloads for their underlying `std::pair` types have been removed.

### Numerical kernel

- Templated the elementary, gamma, beta, Bessel, zeta, Lambert W and Marcum functions on their
  floating-point type. Existing mixed integer/`double` calls remain source-compatible, while
  `long double` now stays in a separate typed path covered by closed-form and recurrence tests.
- Templated factorial and combinatorial return types; calls without an explicit type retain their
  previous `long double` result.

## [1.2.0]

A stabilization release: no API is removed and the file layout is unchanged, but
**sampled values change for the same seed** (see B1–B6, B12, B26 below), and many
probability functions now return different, correct values. Code that stored seeds to
reproduce earlier runs will get different numbers.

The new test suite compares the library with 40-digit [mpmath](https://mpmath.org)
reference values, with the values recorded from v1.1, and with its own internal
coherence and sampling statistics. Every fix below has a regression test in
`tests/unit/test_regressions.cpp` under the same ID.

### Build and portability

- CMake build with `RandLib::RandLib` target, install rules and
  `find_package(RandLib)` support; presets `dev`, `release`, `asan`, `ubsan`.
  `RandLib.pro` (qmake) is deprecated and will be removed in 2.0.
- Builds without warnings on GCC 11–14, Clang 15–18, AppleClang and MSVC 2022 (#18); MSVC's
  narrowing-conversion warnings, which GCC and Clang do not emit, are disabled.
  Duplicate explicit instantiations that MSVC 14.44 diagnoses as C4661 are removed.
  Reference tests document the double-only MSVC incomplete-beta tolerance separately.
  Ziggurat tables and a few constants no longer use `constexpr` initializers that call
  non-`constexpr` math functions. Test constants also avoid names that can trigger MSVC's
  template-instantiation shadow warnings.
- Discrete sampling checks are discovered as one bounded CTest case per distribution. This keeps
  a platform-specific slow or nonterminating sampler from hiding inside a single six-hour test job.
- The library no longer requires the optional C++17 special math functions, which libc++
  does not provide: Bessel I and K, beta and Riemann zeta have portable implementations,
  used automatically when the standard library lacks them
  (`RANDLIB_PORTABLE_SPECIAL_FUNCTIONS=ON` forces them).
- Constants moved to `randlib::consts`; the `M_*` names remain available when the C library
  does not define them. Static builds define `RANDLIB_STATIC` so that Windows consumers do not
  see `dllimport`.
- `RANDLIB_GENERATOR` (CMake) or `RANDLIB_JLKISS64` / `RANDLIB_PCG32` (preprocessor) select the
  engine behind `RandGenerator`; `JLKISS64RAND` is still accepted.
- Repository hygiene: removed `RandLib.pro.user`, LaTeX build files and the dead Travis
  configuration; `NoncentralTRand` moved to `attic/` (it never compiled).

### Random number generation

- **B1** Normal samples no longer depend on other instances: the ziggurat kept a leftover
  exponential in a `static thread_local` variable.
- **B2** JKiss and JLKiss64 seeding expands the seed with SplitMix64 into the whole state;
  seeds that differed only above bit 32 gave identical streams, and the streams of nearby
  seeds were strongly correlated.
- **B3** PCG seeding follows the reference `pcg32_srandom_r`; `PCGRandEngine::Reseed(seed, sequence)`
  selects the stream.
- **B4** Default seeds combine `std::random_device`, a high-resolution clock, the thread id and a
  per-process counter.
- **B5** The 64-bit engine switch works: one uniform generator derives its bit count from the engine.
- **B6** Uniform variates use the full mantissa (52 random bits for `double`, 63 for x87
  `long double`) instead of 32 bits, so inverse-transform tails reach the full range; `float`
  uniforms can no longer round to 1.
- **B9** The normal sign comes from an explicit bit.
- **B12** `BernoulliRand::StandardVariate` no longer caches bits in static variables shared by all
  generators and threads.
- **B26** Composite distributions reseed every generator they use (Laplace and geometric stable
  families, Planck), and member generators get seeds derived with `RandEngine::DeriveSeed`
  instead of `seed + k`, which collided with the seed of a neighbouring object.

### Correctness

- **B7** Stable laws with exponent within 0.01 of 1 and non-zero skewness can be constructed and
  sampled. Their density and distribution function throw the new `randlib::not_implemented`
  until a dedicated algorithm exists (#1).
- **B8** The distribution function of totally skewed stable laws no longer has a floor of
  3.8·10⁻⁵ (the adaptive Simpson integrator missed narrow peaks).
- **B11** Generic quantiles are deterministic and no longer consume the instance's random stream:
  the initial guess came from sampling into a shared static vector (a data race between threads).
  Quantiles are now found by bracketing and the Illinois method to full double precision, in the
  tail with the smaller probability; the old Newton iteration stopped at tolerances of 1e-6–1e-10.
- **B13** `M_LNPI` (log π) was wrong in the 11th digit, which affected every density that uses it.
- **B14** `NormalRand::ThirdMoment` and `FourthMoment` returned wrong values.
- **B15** Classes that could not be used: `BaldingNicholsRand::F` was hidden by a data member,
  `CenteredAsymmetricLaplaceRand` was abstract, `WrappedExponentialRand` and
  `ShiftedAsymmetricLaplaceDistribution` were never instantiated, and the abstract bases had
  protected destructors, so distributions could not be owned through a base pointer.
- **B16** `WrappedExponentialRand`: circular variance, characteristic function at non-integer `t`
  and rate validation. `VonMisesRand` characteristic function at non-integer `t`.
- **B17** Regularized incomplete gamma functions are accurate to about 1e-13 everywhere; they were
  off by up to 1e-3 for large shapes (Gamma, chi-squared, Poisson, Erlang).
- **B18** `logBeta` no longer underflows to −∞ for large arguments (e.g. Beta(600, 600)).
- **B19** Triangular distributions accept a mode equal to a limit.
- **B20** Modified Bessel function of negative and half-integer order: the sinh/cosh cases were
  swapped and the reflection had the wrong sign, which broke the noncentral chi-squared density for
  degrees below 2 and equal to 3.
- **B21** Stable laws: the characteristic function had the skewness term with the wrong sign for
  α ≠ 1, and for α = 1 the density, distribution function and sampler used the opposite sign of the
  S1 location shift (2/π)βγ log γ.
- **B22** Stable densities and distribution functions in the tails: the integrals are split around
  the peak of the integrand, the heavy tail of totally skewed laws with α > 1 no longer returns 0,
  and survival functions are computed directly instead of as 1 − F.
- **B23** Survival functions and tails computed without cancellation or overflow: Cauchy, Student t
  (also near the median), beta prime, Fisher F, inverse Gaussian, logistic, Wigner semicircle.
- **B24** Discrete distributions: logarithmic distribution function, binomial median, categorical
  quantiles, negative binomial survival below the support, Bernoulli log-probability outside the
  support.
- **B25** Invalid parameters throw `std::invalid_argument` instead of being replaced silently
  (log-normal variance, Balding-Nichols, negative binomial), and NaN parameters and non-finite
  locations are rejected throughout.
- **B27** Characteristic functions of heavy-tailed discrete laws (Zeta, Yule) terminate; their
  summation used to need about 10¹⁰ terms. The tail beyond the summed terms is estimated with
  Euler's transformation, so a value takes milliseconds and is accurate to about 1e-16.
- **B28** The Cantor characteristic function ignored `t`.
- **B29** The geometric stable density and distribution function for α = 1 missed the S1 location
  shift of the mixing representation.
- **B30** The Marchenko–Pastur distribution function next to the lower edge of the support.
- **B31** `RademacherRand::Mode` returned a random value; it is now −1, like the median.
- **B32** Poisson and binomial probabilities use Loader's saddle-point form: the sum of logarithms
  lost up to four digits for large arguments, which depended on the platform's `lgamma`.
- **B33** The raised cosine (and Raab–Green) characteristic function at st = π evaluated 0 / 0;
  it is now computed as π² sinc(π − st) / (st (π + st)).
- **B34** Discrete quantiles at a probability level that F reaches exactly no longer depend on the
  last bits of the platform's special functions.
- **B35** Zeta distribution and survival functions use an Euler–Maclaurin tail for large
  arguments instead of summing billions of terms.
- **B36** Continuous and discrete uniform likelihoods negate the sample-count exponent after
  converting it to a signed floating-point value; negating `size_t` previously wrapped and
  produced an enormous positive exponent.
- **B37** Double factorial now evaluates the defining product n(n − 2)…; the previous factorial
  identities shifted by the full input and returned incorrect values for every n > 1.
- Removed the unused private `POISSON` binomial generator selector. Binomial sampling continues to
  use its exact Bernoulli-sum, waiting-time and rejection algorithms rather than an approximation.
- Noncentral chi-squared mode for 2 degrees of freedom; `DegenerateRand::Fit` accepted any sample.
- Numerical integration uses adaptive Gauss–Kronrod (21-point) quadrature with a global error
  estimate instead of adaptive Simpson (#10).
- Quantile searches of stable and geometric stable laws start the far tails from their power-law
  asymptote.

### Known limitations

Tracked by the tolerance overrides in `tests/data/reference_tolerances.csv`:

- characteristic functions of heavy-tailed laws computed by numerical integration are accurate to
  about 1e-3 (#10), the log-normal one to about 1e-4 for large variance (#6);
- noncentral chi-squared upper tails come from 1 − Marcum P (#4);
- stable densities near α = 1 with β ≠ 0 are not implemented (#1);
- geometric stable densities and distribution functions take about 0.1 s per evaluation (up to a
  few seconds in the far tails for α = 1);
- high-precision reference values for the skewed stable and geometric stable laws with α = 1 are
  still being computed; these cases are covered by spot checks against mpmath in the regression
  tests.

## [1.1.0]

Last release before the modernization; tagged in the repository as `v1.1.0-baseline` (6536dad).

# Changelog

All notable changes to RandLib are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [1.2.1] - Unreleased

### Fitting and numerical methods

- **B39** Negative-binomial maximum-likelihood fitting optimizes the positive shape parameter in
  log space, preventing Newton steps from crossing a special-function pole. The Newton helper now
  rejects non-finite states after bounded backtracking instead of performing up to 100,000 futile
  iterations. A three-observation regression uses the 50-digit mpmath score root
  5.2122598849254577174.
- **B40** Rayleigh UMVU scale fitting consistently evaluates the gamma-ratio correction in the log
  domain. The small-sample branch previously passed a logarithm as the scale, while the large-sample
  branch exponentiated half the scale squared and treated a log-ratio expansion as an ordinary
  correction factor.
- **B41** Asymmetric-Laplace entropy now returns
  `1 + log(scale * (asymmetry + 1 / asymmetry))`; the previous `log1p` expression omitted the unit
  entropy contribution and put the scale term inside an unrelated sum.
- **B42** Noncentral chi-squared distributions now accept zero noncentrality, as required by their
  default constructor and central chi-squared limit. Density, direct tails, sampling and the mode use
  exact gamma-limit formulas instead of singular Bessel/Marcum expressions.
- **B43** Bernoulli entropy is now linkable from downstream translation units and returns zero at
  the accepted deterministic endpoints instead of evaluating the indeterminate floating-point
  products `0 * log(0)`.
- **B44** Binomial construction and cached sampling now handle the accepted deterministic
  probabilities without configuring an invalid geometric helper. Binomial and Bernoulli modes at
  `p = 1` are clamped to the upper support endpoint; the ordinary `floor((n + 1) * p)`
  expression otherwise returned `n + 1`. Bulk sampling no longer rounds probabilities within
  `1e-6` of one to a deterministic law.
- **B45** Discrete-uniform characteristic functions now use a reduced-angle sinc form and return
  exactly one at integer-lattice periods instead of dividing zero by zero.
- **B46** Skellam modes now search adjacent probabilities from the median instead of truncating the
  mean, which need not maximize the PMF. A three-point regression uses 60-digit Bessel references.
- **B47** Static Poisson sampling now shares the exact inversion/rejection algorithms used by
  constructed distributions. This replaces its separate exponential-arrival algorithm through rate
  1000 and, more importantly, the inexact normal approximation above that boundary.
- **B48** Geometric and negative-binomial distributions now accept the mathematically valid
  `p = 1` endpoint as a deterministic zero law across probability functions, moments and sampling.
  The improper `p = 0` geometric static sampler is rejected instead of returning the integer limit.
- **B49** Geometric entropy is now reported in nats, consistently with the rest of the API, and
  uses its finite zero-entropy limit at `p = 1`.
- **B50** Bivariate covariance matrices now convert correlation with marginal standard deviations
  rather than variances. Trinomial correlation uses the corresponding square-root formula and is
  undefined when either marginal is deterministic.

## [1.2.0] - 2026-10-17

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
- **B38** Geometric and negative-binomial sampling saturate results that exceed the result integer
  type instead of performing undefined conversions or additions. Binomial waiting-time sampling
  compares against the remaining trials before addition, preventing overflow and a beta-binomial
  hang.
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

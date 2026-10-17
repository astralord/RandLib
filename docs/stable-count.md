# Stable Count distribution

`StableCountRand<Real>` implements the core Stable Count law for `0 < alpha < 1`. Its constructor
uses `(exponent, scale, location)`, and a nonstandard variate is `location + scale * X` for standard
Stable Count `X`. The support is `(location, infinity)`; the density is zero at the finite endpoint.

```cpp
#include <randlib/randlib.hpp>

#include <random>

std::mt19937_64 engine(42);
randlib::StableCountRand<double> count(0.8, 2.0, -1.0);
double draw = count(engine);
double density = count.f(1.0);
double upper_tail = count.S(3.0);
```

This is initially a Tier 2 distribution. The generalized four-parameter family
`(alpha, sigma, d, p)` is not part of this API.

## Definition and stable convention

The definition follows Lihn, *Generalization of the Alpha-Stable Distribution with the Degree of
Freedom*, arXiv:2405.04693v2 (2025), sections 1 and 5.7:

```text
N_alpha(x) = L_alpha(1/x) / (Gamma(1/alpha + 1) x),  x > 0,
E[exp(-s L_alpha)] = exp(-s^alpha).
```

The scale conversion is essential. For `alpha != 1`, RandLib's S1 characteristic function is

```text
exp(-gamma^alpha |t|^alpha (1 - i beta sign(t) tan(pi alpha / 2)) + i mu t).
```

With `beta=1`, choosing

```text
gamma = cos(pi alpha / 2)^(1/alpha),  mu = 0
```

makes the exponent `-(-it)^alpha`; analytic continuation to `t=is` gives the required Laplace
transform `exp(-s^alpha)`. The implementation therefore uses
`StableRand(alpha, 1, gamma, 0, S1)`, not unit stable scale. A focused regression compares both
mappings and fails for the unit-scale convention.

The Wright representation used by the independent oracle is

```text
N_alpha(x) = W_{-alpha,0}(-x^alpha) / Gamma(1/alpha + 1).
```

At `alpha=1/2`, this reduces exactly to `Gamma(shape=3/2, scale=4)`. RandLib uses that identity for
density, probabilities, and sampling. As `alpha` tends to one from below, every standard raw moment
tends to one and the law converges weakly to a point mass at one. Consequently the nonstandard limit
is a point mass at `location + scale`. The continuous endpoint `alpha=1` is rejected.

## Numerical methods

- Density uses the Wright series when cancellation remains bounded, then the reciprocal one-sided
  stable definition. The series uses reflection of `Gamma(-alpha k)`, compensated summation, and a
  finite iteration guard.
- CDF uses the integrated Wright series where safe. Fallback probabilities use the exact change of
  variable
  `F(x) = integral_[1/x,infinity] L_alpha(t)/t dt / Gamma(1/alpha+1)`.
- Small survival probabilities are evaluated directly as
  `integral_[0,1/x] L_alpha(t)/t dt / Gamma(1/alpha+1)`, rather than by subtracting a rounded CDF
  from one. Generic TOMS 748 inversion supplies both quantiles.
- Sampling is correctness-first inverse-CDF sampling from one standard-URBG uniform. The half case
  uses the exact Gamma decomposition. There is no renewal cutoff, discarded large draw, mutable
  static state, or library-owned engine.
- Raw moments use
  `E[X^q] = Gamma((q+1)/alpha)/(Gamma(q+1) Gamma(1/alpha))` for `q > -1`, evaluated through
  log-gamma differences. Mean, variance, skewness, and excess kurtosis include the requested
  location-scale transformation. No unproved closed-form mode is exposed.

The Wright series becomes strongly cancelling in parts of the near-one upper tail. Those regions
deliberately take the slower stable-density quadrature path. Inverse-CDF sampling is correspondingly
slower than specialized transformed samplers, but it has no finite-cutoff bias.

## Independent reference data

`tools/refdata/distributions.py` implements the published Wright series independently of the C++
stable kernel. It uses mpmath 1.3.0 with 40 decimal digits, 120 guard digits, compensated arbitrary-
precision summation, a 20,000-term bound, and precision-doubling confirmation for tiny values. The
half case is evaluated independently with the incomplete gamma function. The committed table covers
`alpha = 0.2, 0.5, 0.8, 0.95`, a location-scale case, both probability tails, quantiles, and moments.

Generation is per-case cached, bounded, resumable, and exactly verifiable:

```sh
uv run --with mpmath==1.3.0 python tools/refdata/generate.py \
  --only StableCountRand --jobs 5 --cache /tmp/randlib-stable-count-refcache
uv run --with mpmath==1.3.0 python tools/refdata/generate.py \
  --only StableCountRand --jobs 5 --cache /tmp/randlib-stable-count-refcache --verify
```

## Fitting

`fit<StableCountRand<double>>(sample)` performs constrained maximum likelihood over exponent, scale,
and location. `FitOptions::fixed` accepts `"exponent"`, `"scale"`, and `"location"`. Initialization
uses exponent `0.5`, a location just below the sample minimum, and moment-scaled spread. Empty,
non-finite, or constant samples fail cleanly.

## Conjugacy research outcome

No `fit_bayes<StableRand>` overload is provided. arXiv:2405.04693v2 twice calls Stable Count a
"conjugate prior", but supplies no normalized likelihood, identifies no random stable parameter,
and gives no posterior or posterior-predictive update. Its proved relations are renewal/counting,
inverse-variable, product, ratio, and subordination identities; those do not establish Bayesian
posterior closure.

The exact half-stable case gives a concrete closure check. If observations have the natural stable-
subordinator likelihood

```text
E[exp(-s Y) | lambda] = exp(-lambda sqrt(s)),
f(y | lambda) = lambda y^(-3/2) exp(-lambda^2/(4y)) / (2 sqrt(pi)),
```

and `lambda` has the standard Stable Count prior, then

```text
pi(lambda) is Gamma(shape=3/2, rate=1/4),
pi(lambda | y_1,...,y_n) is proportional to
  lambda^(n+1/2) exp(-lambda^2 sum(1/y_i)/4 - lambda/4).
```

The quadratic term and changed power show that the posterior is not a location-scale Stable Count
law. Assigning the prior directly to `c=lambda^2` instead gives a Gamma likelihood update whose shape
changes with `n`, again leaving the fixed-shape core Stable Count family. Reciprocal-scale choices
add the corresponding Jacobian and do not repair closure. The paper's random count is a latent
renewal index, not a normalized likelihood parameter for an observed stable sample.

Therefore Phase 6.2 is complete as **investigated but not established**. A future research issue may
define a different likelihood and study closure in the generalized `(alpha, sigma, d, p)` family,
but that requires a new mathematical specification, normalized posterior predictive construction,
and independent quadrature evidence before any public Bayesian API is justified.

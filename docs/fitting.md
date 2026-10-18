# Fitting distributions

RandLib's fitting API is non-mutating: a fit returns `FitResult<Distribution>` containing the fitted
distribution, likelihood and information criteria, uncertainty estimates, convergence metadata, and
the method used.

## Maximum likelihood

```cpp
#include <randlib/randlib.hpp>

#include <vector>

std::vector<double> observations = /* ... */;
auto result = randlib::fit<randlib::GammaRand<double>>(observations);

const auto &model = result.distribution;
double shape = model.GetShape();
double rate = model.GetRate();
double standard_error = result.std_errors.at(0);
```

Pass a `FitOptions` value to control convergence or hold parameters fixed. Parameter names are the
strings returned in `result.parameter_names`; unknown names and impossible values raise `fit_error`.

```cpp
randlib::FitOptions options;
options.fixed = {{"location", 0.0}};
options.tolerance = 1e-9;
auto fitted = randlib::fit<randlib::LevyRand<double>>(observations, options);
```

Closed-form estimators are used where available. Other models use constrained optimization with
positive and bounded parameters transformed to an unconstrained scale. The reported covariance is
the inverse observed information from a numerical Hessian. A fixed or non-identifiable parameter can
have a `NaN` standard error; always check `converged` before using an optimized result.

### Stable laws

Stable fitting returns an S0-parameterized `StableRand`, so the estimated location remains
continuous at `alpha=1`. The fast entry point is McCulloch's quantile estimator:

```cpp
using Stable = randlib::StableRand<double>;
auto estimate = randlib::fit<Stable>(observations, randlib::Method::Quantile);
```

It uses the 5th, 25th, 50th, 75th, and 95th sample percentiles and the published McCulloch lookup
tables for exponent and skewness; S0 scale and location are then matched to theoretical quartiles.
The tables cover exponents down to 0.5. For smaller exponents, use the result as the starting point
for MLE rather than as a final estimate.

`Method::MLE` starts from the same robust estimate and jointly optimizes exponent, skewness, scale,
and location. Stable likelihoods are materially more expensive than ordinary light-tailed models.
Fix known parameters through `FitOptions::fixed`, and use the quantile fit when latency matters.

## Other estimators

`fit_moments<Distribution>(sample)` selects the named method-of-moments entry point. Conjugate
updates return the posterior distribution directly:

```cpp
randlib::NormalInverseGammaRand<double> prior(0.0, 1.0, 2.0, 2.0);
auto posterior = randlib::fit_bayes(observations, prior);
```

Conditional conjugate overloads accept the model's fixed parameter as their final argument. See
`include/randlib/fit/bayes.hpp` for supported prior/model pairs.

## Model selection and checks

```cpp
auto best = randlib::fit_best<randlib::ExponentialRand<double>, randlib::GammaRand<double>,
                              randlib::LogNormalRand<double>>(observations, randlib::Criterion::BIC);
```

`goodness_of_fit(model, sample)` reports Kolmogorov-Smirnov, Anderson-Darling, and
Cramer-von Mises statistics for continuous models. Discrete models use a chi-squared test with sparse
bins combined and fitted parameters accounted for. Information criteria are comparative measures;
they do not replace residual checks or domain knowledge.

## Censoring, truncation, and extremes

Use `fit_censored` with exact, left-, right-, or interval-censored observations and `fit_truncated`
when the observation process excludes values outside known bounds. Both use the same likelihood and
uncertainty pipeline as ordinary MLE.

For heavy tails, `hill_estimator`, `pickands_estimator`, and `hill_plot` provide tail-index
diagnostics. `fit_peaks_over_threshold` fits a generalized Pareto distribution to exceedances above a
chosen threshold. Threshold selection remains a modeling decision: inspect stability over a range of
orders or thresholds rather than accepting one automatic value.

The mutating member functions named `Fit*` are deprecated and will be removed in v3. Define
`RANDLIB_DISABLE_DEPRECATED_FIT_WARNINGS` only as a temporary migration aid.

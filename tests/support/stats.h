#ifndef RANDLIB_TESTS_STATS_H
#define RANDLIB_TESTS_STATS_H

#include <functional>
#include <vector>

/// Statistics used by the tests. They are implemented here, independently of RandLib,
/// so that a bug in the library cannot hide itself in the test.
namespace randlib_test
{

/// Asymptotic p-value of the one-sample Kolmogorov-Smirnov statistic D for n samples,
/// using Stephens' small-sample correction.
double kolmogorovSmirnovPValue(double D, size_t n);

/// Kolmogorov-Smirnov statistic of a sorted sample against a cdf.
double kolmogorovSmirnovStatistic(const std::vector<double> &sorted, const std::function<double(double)> &cdf);

/// Upper tail P(X > x) of the chi-squared distribution with k degrees of freedom.
double chiSquaredSurvival(double x, double k);

/// Integral of f over [a, b] by adaptive Gauss-Kronrod (7, 15) quadrature.
double integrate(const std::function<double(double)> &f, double a, double b, double tolerance = 1e-12);

/// True when |actual - expected| <= tolerance * max(|actual|, |expected|), or both are equal
/// (including equal infinities and both NaN).
bool relativelyClose(double actual, double expected, double tolerance);

/// Relative difference |actual - expected| / |expected| (absolute when expected is 0).
double relativeError(double actual, double expected);

/// Runs body(0), ..., body(count - 1) on a pool of worker threads.
/// Catch2 assertions are not thread-safe: collect results and assert afterwards.
void parallelFor(size_t count, const std::function<void(size_t)> &body);

}

#endif // RANDLIB_TESTS_STATS_H

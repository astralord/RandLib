// Quadrature, root finding and optimization in randlib::math.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <numbers>
#include <span>
#include <vector>

using namespace randlib;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Numerical functions preserve long-double arithmetic", "[numerics]")
{
    using Real = long double;
    auto checkClose = [](Real actual, Real expected, Real relativeTolerance = 2e-16L) {
        CHECK(std::fabs(actual - expected) <= relativeTolerance * std::max(Real{1}, std::fabs(expected)));
    };

    STATIC_CHECK(std::same_as<decltype(math::softplus(Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::logBeta(Real{1}, Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::pgamma(Real{1}, Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::besselI(Real{0}, Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::factorial(Real{1})), Real>);
    STATIC_CHECK(std::same_as<decltype(math::factorial(1)), long double>);
    STATIC_CHECK(std::same_as<decltype(math::doubleFactorial<Real>(1)), Real>);
    STATIC_CHECK(std::same_as<decltype(math::lfact<Real>(1)), Real>);

    checkClose(math::atan(Real{2}), std::atan(Real{2}));
    checkClose(math::softplus(Real{-2}), std::log1p(std::exp(Real{-2})));
    checkClose(math::log1mexp(Real{-2}), std::log1p(-std::exp(Real{-2})));
    checkClose(math::logexpm1l(Real{2}), std::log(std::expm1(Real{2})));
    checkClose(math::log2mexp(Real{-2}), std::log(Real{2} - std::exp(Real{-2})));

    constexpr Real inverseArgument = 0.375L;
    checkClose(math::erfinv(std::erf(inverseArgument)), inverseArgument, 2e-15L);
    checkClose(math::erfcinv(std::erfc(inverseArgument)), inverseArgument, 2e-15L);
    checkClose(math::erfcx(Real{2}), std::exp(Real{4}) * std::erfc(Real{2}), 2e-15L);
    checkClose(math::harmonicNumber(Real{1}, 4), Real{25} / 12, 1e-14L);

    checkClose(math::factorial(Real{10}), Real{3628800});
    checkClose(math::lfact<Real>(10), std::lgamma(Real{11}));
    checkClose(math::ldfact<Real>(9), std::log(Real{945}));
    checkClose(math::binom<Real>(20, 6), Real{38760}, 2e-15L);

    constexpr Real x = 2.25L;
    checkClose(math::digamma(x + 1), math::digamma(x) + 1 / x);
    checkClose(math::trigamma(x + 1), math::trigamma(x) - 1 / (x * x));
    checkClose(math::ibeta(Real{0.25}, Real{2}, Real{3}), Real{0.26171875}, 2e-15L);
    checkClose(math::pgamma(Real{1}, x), -std::expm1(-x), 2e-15L);
    checkClose(math::pgamma(Real{2}, x) + math::qgamma(Real{2}, x), Real{1}, 2e-15L);

    constexpr Real besselArgument = 1.25L;
    Real expectedBesselI = std::sqrt(Real{2} / (std::numbers::pi_v<Real> * besselArgument)) * std::sinh(besselArgument);
    checkClose(math::portable::besselI(Real{0.5}, besselArgument), expectedBesselI, 2e-15L);
    checkClose(math::portable::riemannZeta(Real{2}), std::numbers::pi_v<Real> * std::numbers::pi_v<Real> / 6, 2e-15L);

    Real w = math::W0Lambert(Real{3});
    checkClose(w * std::exp(w), Real{3}, 2e-15L);
    Real marcumP = math::MarcumP(Real{1}, Real{2}, Real{3});
    Real marcumQ = math::MarcumQ(Real{1}, Real{2}, Real{3});
    CHECK(marcumP >= 0);
    CHECK(marcumP <= 1);
    checkClose(marcumP + marcumQ, Real{1});
}

TEST_CASE("Adaptive quadrature on finite and infinite intervals", "[numerics]")
{
    CHECK_THAT(math::integrate([](double x) { return std::exp(x); }, 0, 1), WithinRel(std::numbers::e - 1, 1e-13));
    CHECK_THAT(math::integrate([](double x) { return 1 / std::sqrt(x); }, 0, 1, 1e-12, 4096), WithinRel(2.0, 1e-10));
    CHECK_THAT(math::integrate([](double x) { return std::exp(-x * x); }, -INFINITY, INFINITY),
               WithinRel(std::sqrt(std::numbers::pi), 1e-12));
    CHECK_THAT(math::integrate([](double x) { return std::exp(-x); }, 2, INFINITY), WithinRel(std::exp(-2.0), 1e-12));
    CHECK_THAT(math::integrate([](double x) { return 1 / (1 + x * x); }, -INFINITY, 0),
               WithinRel(std::numbers::pi / 2, 1e-12));
    CHECK(math::integrate([](double x) { return x; }, 1, 0) == -math::integrate([](double x) { return x; }, 0, 1));
    CHECK_THAT(static_cast<double>(math::integrate<long double>([](long double x) { return x * x; }, 0, 3)),
               WithinRel(9.0, 1e-15));
}

TEST_CASE("TOMS 748 finds bracketed roots to full precision", "[numerics]")
{
    double root = math::findRoot([](double x) { return std::cos(x) - x; }, 0, 1);
    CHECK_THAT(root, WithinRel(0.73908513321516064166, 4e-16));
    double cubic = math::findRoot([](double x) { return x * x * x - 2; }, 0, 10);
    CHECK_THAT(cubic, WithinRel(std::cbrt(2.0), 4e-16));
    /// a steep function with a root near one end
    double steep = math::findRoot([](double x) { return std::exp(20 * x) - 1e-3; }, -5, 5);
    CHECK_THAT(steep, WithinRel(std::log(1e-3) / 20, 1e-14));
    /// a discontinuous function: the bracket closes on the jump
    double jump = math::findRoot([](double x) { return x < 0.3 ? -1.0 : 1.0; }, 0, 1);
    CHECK_THAT(jump, WithinAbs(0.3, 1e-15));
    CHECK_THROWS_AS(math::findRoot([](double x) { return x * x + 1; }, -1, 1), std::invalid_argument);
}

TEST_CASE("Multivariate minimization", "[numerics]")
{
    auto rosenbrock = [](std::span<const double> x) {
        return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
    };
    auto simplex = math::nelderMead(rosenbrock, {-1.2, 1.0});
    CHECK(simplex.converged);
    CHECK_THAT(simplex.x[0], WithinAbs(1.0, 1e-6));
    CHECK_THAT(simplex.x[1], WithinAbs(1.0, 1e-6));

    auto quasiNewton = math::bfgs(rosenbrock, {-1.2, 1.0});
    CHECK(quasiNewton.converged);
    CHECK_THAT(quasiNewton.x[0], WithinAbs(1.0, 1e-5));
    CHECK_THAT(quasiNewton.x[1], WithinAbs(1.0, 1e-5));

    auto quadratic = [](std::span<const double> x) {
        double sum = 0;
        for (size_t i = 0; i < x.size(); ++i)
            sum += (i + 1) * (x[i] - i) * (x[i] - i);
        return sum;
    };
    auto result = math::bfgs(quadratic, std::vector<double>(5, 10.0));
    for (size_t i = 0; i < 5; ++i)
        CHECK_THAT(result.x[i], WithinAbs(static_cast<double>(i), 1e-6));
}

TEST_CASE("Numerical Hessian recovers quadratic curvature", "[numerics]")
{
    auto quadratic = [](std::span<const double> x) {
        return 2 * x[0] * x[0] + 3 * x[0] * x[1] + 5 * x[1] * x[1] + 7 * x[0] - 11 * x[1];
    };
    auto hessian = math::numericalHessian(quadratic, std::vector<double>{1.25, -0.75});
    REQUIRE(hessian.size() == 2);
    REQUIRE(hessian[0].size() == 2);
    CHECK_THAT(hessian[0][0], WithinAbs(4.0, 1e-6));
    CHECK_THAT(hessian[0][1], WithinAbs(3.0, 1e-6));
    CHECK_THAT(hessian[1][0], WithinAbs(3.0, 1e-6));
    CHECK_THAT(hessian[1][1], WithinAbs(10.0, 1e-6));
}

TEST_CASE("Exponential-family measures support scalar and pair parameters", "[numerics]")
{
    NormalRand<double> normal(1.5, 2.0);
    CHECK_THAT(normal.LogProbabilityMeasure(0.25), WithinRel(normal.logf(0.25), 1e-15));

    constexpr double otherMean = -0.5;
    constexpr double otherVariance = 3.0;
    const double expectedDivergence =
        0.5 * (std::log(otherVariance / normal.Variance()) +
               (normal.Variance() + std::pow(normal.Mean() - otherMean, 2)) / otherVariance - 1);
    CHECK_THAT(normal.KullbackLeiblerDivergence({otherMean, otherVariance}), WithinRel(expectedDivergence, 1e-15));

    PoissonRand<int> poisson(3.5);
    CHECK_THAT(poisson.LogProbabilityMeasure(2), WithinRel(poisson.logP(2), 1e-15));
}

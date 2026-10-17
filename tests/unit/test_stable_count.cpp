#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <randlib/randlib.hpp>

#include <cmath>
#include <limits>
#include <random>

using Catch::Approx;
using namespace randlib;

TEST_CASE("Stable Count validates parameters and support", "[stable-count]")
{
    CHECK_THROWS_AS(StableCountRand<double>(0), std::invalid_argument);
    CHECK_THROWS_AS(StableCountRand<double>(1), std::invalid_argument);
    CHECK_THROWS_AS(StableCountRand<double>(NAN), std::invalid_argument);
    CHECK_THROWS_AS(StableCountRand<double>(0.5, 0), std::invalid_argument);
    CHECK_THROWS_AS(StableCountRand<double>(0.5, 1, INFINITY), std::invalid_argument);

    StableCountRand<double> distribution(0.05, 2, -3);
    CHECK(distribution.GetExponent() == 0.05);
    CHECK(distribution.GetScale() == 2);
    CHECK(distribution.GetLocation() == -3);
    CHECK(distribution.SupportType() == Support::RightSemiFinite);
    CHECK(distribution.MinValue() == -3);
    CHECK(std::isinf(distribution.MaxValue()));
    CHECK(distribution.f(-3) == 0);
    CHECK(distribution.logf(-3) == -INFINITY);
    CHECK(distribution.F(-3) == 0);
    CHECK(distribution.S(-3) == 1);
    CHECK(distribution.Quantile(0) == -3);
    CHECK(std::isinf(distribution.Quantile(1)));
    CHECK(distribution.f(-1) == Approx(3.7864403055304811e-21).epsilon(2e-13));
    CHECK(distribution.F(-1) == Approx(7.54400211146346e-21).epsilon(2e-13));
}

TEST_CASE("Stable Count uses the positive-stable S1 scale convention", "[stable-count]")
{
    constexpr double alpha = 0.8;
    constexpr double x = 1.2;
    StableCountRand<double> distribution(alpha);
    const double stableScale = std::pow(std::cos(randlib::consts::half_pi * alpha), 1 / alpha);
    StableRand<double> positiveStable(alpha, 1, stableScale, 0, StableRand<double>::Parameterization::S1);
    const double mapped = positiveStable.f(1 / x) / (std::tgamma(1 / alpha + 1) * x);
    CHECK(distribution.f(x) == Approx(mapped).epsilon(2e-12));

    StableRand<double> wrongScale(alpha, 1, 1, 0, StableRand<double>::Parameterization::S1);
    const double wronglyMapped = wrongScale.f(1 / x) / (std::tgamma(1 / alpha + 1) * x);
    CHECK(std::fabs(wronglyMapped / distribution.f(x) - 1) > 0.1);
}

TEST_CASE("Half Stable Count is exactly Gamma three-halves with scale four", "[stable-count]")
{
    StableCountRand<double> stableCount;
    const auto gamma = GammaRand<double>::FromShapeScale(1.5, 4);
    for (double x : {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0}) {
        CHECK(stableCount.logf(x) == Approx(gamma.logf(x)).margin(2e-14));
        CHECK(stableCount.f(x) == Approx(gamma.f(x)).epsilon(2e-14));
        CHECK(stableCount.F(x) == Approx(gamma.F(x)).epsilon(2e-14));
        CHECK(stableCount.S(x) == Approx(gamma.S(x)).epsilon(2e-14));
    }
    for (double p : {1e-8, 0.01, 0.25, 0.5, 0.9, 1 - 1e-8})
        CHECK(stableCount.Quantile(p) == Approx(gamma.Quantile(p)).epsilon(2e-10));

    std::mt19937_64 left(1234), right(1234);
    for (int i = 0; i < 100; ++i)
        CHECK(stableCount(left) == gamma(right));
}

TEST_CASE("Stable Count moments transform under location and scale", "[stable-count]")
{
    StableCountRand<double> standard;
    CHECK(standard.StandardMoment(0) == Approx(1).epsilon(2e-15));
    CHECK(standard.Mean() == Approx(6).epsilon(2e-15));
    CHECK(standard.SecondMoment() == Approx(60).epsilon(2e-15));
    CHECK(standard.Variance() == Approx(24).epsilon(2e-15));
    CHECK(standard.Skewness() == Approx(2 / std::sqrt(1.5)).epsilon(2e-14));
    CHECK(standard.ExcessKurtosis() == Approx(4).epsilon(2e-14));

    StableCountRand<double> transformed(0.5, 3, -2);
    CHECK(transformed.Mean() == Approx(16).epsilon(2e-15));
    CHECK(transformed.Variance() == Approx(216).epsilon(2e-15));
    CHECK(transformed.SecondMoment() == Approx(472).epsilon(2e-15));
    CHECK(transformed.Skewness() == Approx(standard.Skewness()).epsilon(2e-15));
    CHECK(transformed.ExcessKurtosis() == Approx(standard.ExcessKurtosis()).epsilon(2e-15));
}

TEST_CASE("Stable Count probability functions and quantiles preserve both tails", "[stable-count]")
{
    StableCountRand<double> distribution(0.8);
    for (double x : {0.1, 0.5, 1.0, 1.2, 2.0}) {
        CHECK(distribution.f(x) == Approx(std::exp(distribution.logf(x))).epsilon(2e-14));
        CHECK(distribution.F(x) + distribution.S(x) == Approx(1).margin(2e-14));
    }
    CHECK(distribution.F(1) == Approx(0.202918804027770297957828129476786).epsilon(2e-13));
    CHECK(distribution.S(5) > 0);
    CHECK(distribution.F(5) == 1);
    for (double p : {0.01, 0.25, 0.5, 0.75, 0.99}) {
        const double quantile = distribution.Quantile(p);
        CHECK(distribution.F(quantile) == Approx(p).margin(2e-8));
    }

    const double normalization =
        math::integrate<double>([&distribution](double x) { return distribution.f(x); }, 0, INFINITY, 2e-9, 512);
    CHECK(normalization == Approx(1).margin(2e-8));
}

TEST_CASE("Stable Count evaluates moderate upper tails directly", "[stable-count]")
{
    StableCountRand<double> distribution(0.2);
    constexpr double upperTailQuantile = 171301.11941043584;
    CHECK(distribution.S(upperTailQuantile) == Approx(0.01).margin(2e-14));
    CHECK(distribution.F(upperTailQuantile) + distribution.S(upperTailQuantile) == Approx(1).margin(2e-14));
}

TEST_CASE("Stable Count approaches its degenerate alpha-one limit", "[stable-count]")
{
    StableCountRand<double> distribution(0.999);
    CHECK(distribution.Mean() == Approx(1).margin(0.01));
    CHECK(distribution.Variance() < 0.01);
    CHECK(distribution.StandardMoment(3) == Approx(1).margin(0.05));
}

TEST_CASE("Stable Count sampling accepts standard engines and has no shared state", "[stable-count]")
{
    StableCountRand<double> distribution(0.8, 2, -1);
    std::minstd_rand first(42), second(42), unrelated(7);
    for (int i = 0; i < 8; ++i) {
        static_cast<void>(distribution(unrelated));
        CHECK(distribution(first) == distribution(second));
    }

    StableCountRand<long double> extended(0.5L, 2, -1);
    CHECK(extended.f(3) == Approx(0.5 * GammaRand<long double>::FromShapeScale(1.5, 4).f(2)).epsilon(2e-14));
}

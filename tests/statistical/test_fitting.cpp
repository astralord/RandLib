// Parameter recovery and likelihood maximality of the maximum-likelihood Fit* methods.
// Each test draws a sample from known parameters (fixed seed), fits it, checks that the
// estimates are close to the truth and that no small perturbation of the estimates
// increases the log-likelihood.

#include <catch2/catch_test_macros.hpp>

#include "RandLib.h"

#include <cmath>
#include <functional>
#include <vector>

namespace
{

constexpr size_t N = 20000;
constexpr double RECOVERY = 0.05;
constexpr double STEP = 1e-3;

template <class D, class T>
std::vector<T> draw(D distribution, unsigned long seed, size_t n = N)
{
    std::vector<T> sample(n);
    distribution.Reseed(seed);
    distribution.Sample(sample);
    return sample;
}

/// Perturbs each free parameter by ±STEP (relative) and requires the log-likelihood not to increase
void requireLocalMaximum(const std::function<double(const std::vector<double> &)> &logLikelihood,
                         const std::vector<double> &estimate, const std::vector<bool> &free)
{
    double best = logLikelihood(estimate);
    REQUIRE(std::isfinite(best));
    for (size_t i = 0; i < estimate.size(); ++i) {
        if (!free[i])
            continue;
        for (double sign : {-1.0, 1.0}) {
            std::vector<double> moved = estimate;
            moved[i] = (moved[i] == 0) ? sign * STEP : moved[i] * (1 + sign * STEP);
            double value;
            try {
                value = logLikelihood(moved);
            } catch (const std::invalid_argument &) {
                continue;
            }
            INFO("parameter " << i << " moved from " << estimate[i] << " to " << moved[i]);
            CHECK(value <= best + 1e-9 * std::fabs(best));
        }
    }
}

void requireClose(double estimate, double truth, double tolerance = RECOVERY)
{
    INFO("estimate " << estimate << ", truth " << truth);
    CHECK(std::fabs(estimate - truth) <= tolerance * std::max(std::fabs(truth), 1.0));
}

}

TEST_CASE("Normal maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<NormalRand<double>, double>(NormalRand<double>(1.5, 4), 1);
    NormalRand<double> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetLocation(), 1.5);
    requireClose(static_cast<double>(fitted.Variance()), 4);
    requireLocalMaximum([&](const std::vector<double> &p) { return NormalRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                        {fitted.GetLocation(), static_cast<double>(fitted.Variance())}, {true, true});
}

TEST_CASE("Log-normal maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<LogNormalRand<double>, double>(LogNormalRand<double>(0.5, 0.25), 2);
    LogNormalRand<double> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetLocation(), 0.5);
    requireClose(fitted.GetScale() * fitted.GetScale(), 0.25);
    requireLocalMaximum([&](const std::vector<double> &p) { return LogNormalRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                        {fitted.GetLocation(), fitted.GetScale() * fitted.GetScale()}, {true, true});
}

TEST_CASE("Exponential maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<ExponentialRand<double>, double>(ExponentialRand<double>(2.5), 3);
    ExponentialRand<double> fitted;
    fitted.FitRate(sample);
    requireClose(fitted.GetRate() / 2.5, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return ExponentialRand<double>(p[0]).LogLikelihoodFunction(sample); },
                        {fitted.GetRate()}, {true});
}

TEST_CASE("Gamma maximum-likelihood fit", "[fitting]")
{
    for (auto [shape, rate] : {std::pair{0.3, 2.0}, std::pair{2.5, 0.5}, std::pair{50.0, 3.0}}) {
        INFO("Gamma(" << shape << ", " << rate << ")");
        auto sample = draw<GammaRand<double>, double>(GammaRand<double>(shape, rate), 4);
        GammaRand<double> fitted;
        fitted.Fit(sample);
        requireClose(fitted.GetShape() / shape, 1);
        requireClose(fitted.GetRate() / rate, 1);
        requireLocalMaximum([&](const std::vector<double> &p) { return GammaRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                            {fitted.GetShape(), fitted.GetRate()}, {true, true});
    }
}

TEST_CASE("Beta maximum-likelihood fit", "[fitting]")
{
    for (auto [a, b] : {std::pair{0.5, 0.5}, std::pair{2.0, 5.0}, std::pair{30.0, 20.0}}) {
        INFO("Beta(" << a << ", " << b << ")");
        auto sample = draw<BetaRand<double>, double>(BetaRand<double>(a, b), 5);
        BetaRand<double> fitted;
        fitted.FitShapes(sample);
        requireClose(fitted.GetAlpha() / a, 1);
        requireClose(fitted.GetBeta() / b, 1);
        requireLocalMaximum([&](const std::vector<double> &p) { return BetaRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                            {fitted.GetAlpha(), fitted.GetBeta()}, {true, true});
    }
}

TEST_CASE("Beta-prime maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<BetaPrimeRand<double>, double>(BetaPrimeRand<double>(2, 3), 6);
    BetaPrimeRand<double> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetAlpha() / 2, 1);
    requireClose(fitted.GetBeta() / 3, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return BetaPrimeRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                        {fitted.GetAlpha(), fitted.GetBeta()}, {true, true});
}

TEST_CASE("Pareto maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<ParetoRand<double>, double>(ParetoRand<double>(3, 2), 7);
    ParetoRand<double> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetShape() / 3, 1);
    requireClose(fitted.GetScale() / 2, 1, 1e-3);
    requireLocalMaximum([&](const std::vector<double> &p) { return ParetoRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                        {fitted.GetShape(), fitted.GetScale()}, {true, true});
}

TEST_CASE("Weibull maximum-likelihood fit of the scale", "[fitting]")
{
    auto sample = draw<WeibullRand<double>, double>(WeibullRand<double>(2, 1.5), 8);
    WeibullRand<double> fitted(1, 1.5);
    fitted.FitScale(sample);
    requireClose(fitted.GetScale() / 2, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return WeibullRand<double>(p[0], 1.5).LogLikelihoodFunction(sample); },
                        {fitted.GetScale()}, {true});
}

TEST_CASE("Laplace maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<LaplaceRand<double>, double>(LaplaceRand<double>(2, 0.5), 9);
    LaplaceRand<double> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetShift(), 2);
    requireClose(fitted.GetScale() / 0.5, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return LaplaceRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                        {fitted.GetShift(), fitted.GetScale()}, {true, true});
}

TEST_CASE("Logistic maximum-likelihood fit of the location", "[fitting]")
{
    auto sample = draw<LogisticRand<double>, double>(LogisticRand<double>(3, 0.4), 10);
    LogisticRand<double> fitted(0, 0.4);
    fitted.FitLocation(sample);
    requireClose(fitted.GetLocation(), 3);
    requireLocalMaximum([&](const std::vector<double> &p) { return LogisticRand<double>(p[0], 0.4).LogLikelihoodFunction(sample); },
                        {fitted.GetLocation()}, {true});
}

TEST_CASE("Levy maximum-likelihood fit of the scale", "[fitting]")
{
    auto sample = draw<LevyRand<double>, double>(LevyRand<double>(0, 3), 11);
    LevyRand<double> fitted(0, 1);
    fitted.FitScale(sample);
    requireClose(fitted.GetScale() / 3, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return LevyRand<double>(0, p[0]).LogLikelihoodFunction(sample); },
                        {fitted.GetScale()}, {true});
}

TEST_CASE("Rayleigh maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<RayleighRand<double>, double>(RayleighRand<double>(0.4), 12);
    RayleighRand<double> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetScale() / 0.4, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return RayleighRand<double>(p[0]).LogLikelihoodFunction(sample); },
                        {fitted.GetScale()}, {true});
}

TEST_CASE("Uniform maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<UniformRand<double>, double>(UniformRand<double>(-2, 5), 13);
    UniformRand<double> fitted;
    fitted.Fit(sample);
    requireClose(fitted.MinValue(), -2, 1e-3);
    requireClose(fitted.MaxValue(), 5, 1e-3);
    requireLocalMaximum([&](const std::vector<double> &p) { return UniformRand<double>(p[0], p[1]).LogLikelihoodFunction(sample); },
                        {fitted.MinValue(), fitted.MaxValue()}, {true, true});
}

TEST_CASE("Poisson maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<PoissonRand<int>, int>(PoissonRand<int>(3.5), 14);
    PoissonRand<int> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetRate() / 3.5, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return PoissonRand<int>(p[0]).LogLikelihoodFunction(sample); },
                        {fitted.GetRate()}, {true});
}

TEST_CASE("Binomial maximum-likelihood fit of the probability", "[fitting]")
{
    auto sample = draw<BinomialRand<int>, int>(BinomialRand<int>(10, 0.3), 15);
    BinomialRand<int> fitted(10, 0.5);
    fitted.FitProbability(sample);
    requireClose(fitted.GetProbability() / 0.3, 1);
    requireLocalMaximum([&](const std::vector<double> &p) { return BinomialRand<int>(10, p[0]).LogLikelihoodFunction(sample); },
                        {fitted.GetProbability()}, {true});
}

TEST_CASE("Negative binomial maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<NegativeBinomialRand<int>, int>(NegativeBinomialRand<int>(2.5, 0.4), 16);
    NegativeBinomialRand<int> fitted;
    fitted.Fit(sample);
    requireClose(fitted.GetNumber() / 2.5, 1, 0.1);
    requireClose(fitted.GetProbability() / 0.4, 1, 0.1);
    requireLocalMaximum([&](const std::vector<double> &p) { return NegativeBinomialRand<int>(p[0], p[1]).LogLikelihoodFunction(sample); },
                        {static_cast<double>(fitted.GetNumber()), fitted.GetProbability()}, {true, true});
}

TEST_CASE("Uniform discrete maximum-likelihood fit", "[fitting]")
{
    auto sample = draw<UniformDiscreteRand<int>, int>(UniformDiscreteRand<int>(-5, 5), 17);
    UniformDiscreteRand<int> fitted;
    fitted.Fit(sample);
    CHECK(fitted.MinValue() == -5);
    CHECK(fitted.MaxValue() == 5);
}

TEST_CASE("Degenerate fit", "[fitting]")
{
    std::vector<double> sample(10, 2.5);
    DegenerateRand fitted;
    fitted.Fit(sample);
    CHECK(fitted.Mean() == 2.5);
    sample.push_back(3);
    CHECK_THROWS_AS(fitted.Fit(sample), std::invalid_argument);
}

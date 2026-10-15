#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include <cmath>
#include <vector>

using namespace randlib;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Unified fitting reports likelihood and information criteria", "[fit]")
{
    const std::vector<double> sample = {1, 2, 3, 4};
    auto result = fit<NormalRand<double>>(sample, Method::MLE);

    CHECK(result.converged);
    CHECK(result.iterations == 0);
    CHECK(result.method == Method::MLE);
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(2.5, 1e-14));
    CHECK_THAT(static_cast<double>(result.distribution.Variance()), WithinAbs(1.25, 1e-14));
    CHECK(result.parameter_names == std::vector<std::string>{"location", "variance"});
    REQUIRE(result.parameters.size() == 2);
    CHECK_THAT(result.parameters[0], WithinAbs(2.5, 1e-14));
    CHECK_THAT(result.parameters[1], WithinAbs(1.25, 1e-14));
    CHECK(result.std_errors.size() == 2);
    CHECK(std::isfinite(result.std_errors[0]));
    CHECK(std::isfinite(result.std_errors[1]));
    CHECK_THAT(result.aic, WithinAbs(4 - 2 * result.log_likelihood, 1e-14));
    CHECK_THAT(result.bic, WithinAbs(2 * std::log(4.0) - 2 * result.log_likelihood, 1e-14));
}

TEST_CASE("Unified fitting supports fixed parameters", "[fit]")
{
    const std::vector<double> sample = {1, 2, 3, 4};
    FitOptions options;
    options.fixed = {{"variance", 2.0}};
    auto result = fit<NormalRand<double>>(sample, options);

    CHECK(result.converged);
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(2.5, 1e-7));
    CHECK_THAT(static_cast<double>(result.distribution.Variance()), WithinAbs(2.0, 1e-14));
    CHECK(result.fixed == std::vector<bool>{false, true});
    CHECK(result.std_errors[1] == 0);
    CHECK_THAT(result.aic, WithinAbs(2 - 2 * result.log_likelihood, 1e-14));
}

TEST_CASE("Unified fitting supports discrete distributions", "[fit]")
{
    const std::vector<int> sample = {0, 1, 2, 3};
    auto result = fit<PoissonRand<int>>(sample);

    CHECK(result.converged);
    CHECK_THAT(result.distribution.GetRate(), WithinAbs(1.5, 1e-14));
    CHECK(result.parameter_names == std::vector<std::string>{"rate"});
    CHECK(result.std_errors[0] > 0);
}

TEST_CASE("Unified fitting rejects invalid options", "[fit]")
{
    const std::vector<double> sample = {1, 2, 3};
    FitOptions options;
    options.fixed = {{"unknown", 1.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), std::invalid_argument);
    options.fixed = {{"variance", -1.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), std::invalid_argument);
    options.fixed = {{"variance", 0.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), std::invalid_argument);
    options.method = Method::UMVU;
    options.fixed = {{"variance", 1.0}};
    CHECK_THROWS_AS(fit<NormalRand<double>>(sample, options), fit_error);
    CHECK_THROWS_AS(fit<NormalRand<double>>(std::span<const double>{}), fit_error);
}

TEST_CASE("Unified fitting exposes existing continuous estimators", "[fit]")
{
    const std::vector<double> positive = {0.2, 0.5, 0.8, 1.1, 1.7, 2.4, 3.2};
    const std::vector<double> unit = {0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9};
    const std::vector<double> real = {-3, -1, -0.5, 0, 0.25, 2, 5};

    CHECK(fit<LogNormalRand<double>>(positive).converged);
    CHECK(fit<ExponentialRand<double>>(positive).converged);
    CHECK(fit<GammaRand<double>>(positive).converged);
    CHECK(fit<BetaRand<double>>(unit).converged);
    CHECK(fit<ArcsineRand<double>>(unit).converged);
    CHECK(fit<BetaPrimeRand<double>>(positive).converged);
    CHECK(fit<ParetoRand<double>>(positive).converged);
    CHECK(fit<LaplaceRand<double>>(real).converged);
    CHECK(fit<RayleighRand<double>>(positive).converged);
    CHECK(fit<UniformRand<double>>(real).converged);

    const std::vector<double> constant(5, 2.5);
    auto degenerate = fit<DegenerateRand>(constant);
    CHECK_THAT(degenerate.distribution.GetValue(), WithinAbs(2.5, 0));
}

TEST_CASE("Unified fitting handles distributions with fixed structural parameters", "[fit]")
{
    const std::vector<int> binomialSample = {1, 2, 2, 3, 4, 3, 2, 1};
    FitOptions binomialOptions;
    binomialOptions.fixed = {{"number", 10}};
    auto binomial = fit<BinomialRand<int>>(binomialSample, binomialOptions);
    CHECK(binomial.converged);
    CHECK(binomial.distribution.GetNumber() == 10);
    CHECK_THAT(binomial.distribution.GetProbability(), WithinAbs(0.225, 1e-8));

    const std::vector<double> levySample = {0.1, 0.25, 0.5, 1, 2, 4, 8};
    FitOptions levyOptions;
    levyOptions.fixed = {{"location", 0}};
    auto levy = fit<LevyRand<double>>(levySample, levyOptions);
    CHECK(levy.converged);
    CHECK(levy.distribution.GetLocation() == 0);

    FitOptions invalidBinomial;
    invalidBinomial.fixed = {{"number", 3.5}};
    CHECK_THROWS_AS(fit<BinomialRand<int>>(binomialSample, invalidBinomial), std::invalid_argument);
    CHECK_THROWS_AS(fit<BinomialRand<int>>(binomialSample), fit_error);
}

TEST_CASE("Unified fitting covers discrete closed-form estimators", "[fit]")
{
    const std::vector<int> overdispersed = {0, 0, 1, 1, 2, 3, 6, 12};
    auto negativeBinomial = fit<PolyaRand<int>>(overdispersed);
    CHECK(negativeBinomial.converged);
    CHECK(negativeBinomial.distribution.GetNumber() > 0);
    CHECK(negativeBinomial.distribution.GetProbability() > 0);
    CHECK(negativeBinomial.distribution.GetProbability() < 1);

    const std::vector<int> bounded = {-2, -1, 0, 2, 5};
    auto uniform = fit<UniformDiscreteRand<int>>(bounded);
    CHECK(uniform.converged);
    CHECK(uniform.distribution.MinValue() == -2);
    CHECK(uniform.distribution.MaxValue() == 5);
}

TEST_CASE("Generic MLE recovers Weibull and logistic parameters", "[fit]")
{
    constexpr int sampleSize = 400;
    std::vector<double> weibullSample;
    std::vector<double> logisticSample;
    weibullSample.reserve(sampleSize);
    logisticSample.reserve(sampleSize);
    const WeibullRand<double> weibullReference(2.5, 1.7);
    const LogisticRand<double> logisticReference(-1.25, 0.8);
    for (int i = 0; i < sampleSize; ++i) {
        const double probability = (i + 0.5) / sampleSize;
        weibullSample.push_back(weibullReference.Quantile(probability));
        logisticSample.push_back(logisticReference.Quantile(probability));
    }

    auto weibull = fit<WeibullRand<double>>(weibullSample);
    REQUIRE(weibull.converged);
    CHECK_THAT(weibull.distribution.GetScale(), WithinRel(2.5, 0.01));
    CHECK_THAT(weibull.distribution.GetShape(), WithinRel(1.7, 0.01));

    auto logistic = fit<LogisticRand<double>>(logisticSample);
    REQUIRE(logistic.converged);
    CHECK_THAT(logistic.distribution.GetLocation(), WithinAbs(-1.25, 1e-8));
    CHECK_THAT(logistic.distribution.GetScale(), WithinRel(0.8, 0.01));
}

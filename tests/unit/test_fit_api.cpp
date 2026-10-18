#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <randlib/randlib.hpp>

#include "data.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

using namespace randlib;
using namespace randlib_test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

struct ScipyFitCase
{
    std::string distribution;
    std::vector<double> sample;
    std::map<std::string, double, std::less<>> parameters;
};

std::map<std::string, ScipyFitCase, std::less<>> loadScipyFitCases()
{
    std::map<std::string, ScipyFitCase, std::less<>> cases;
    std::ifstream samples(dataDirectory() + "/fitting/scipy_samples.csv");
    std::string line;
    std::getline(samples, line);
    while (std::getline(samples, line)) {
        auto fields = splitCsvLine(line);
        cases[fields.at(0)].sample.push_back(parseNumber(fields.at(1)));
    }

    std::ifstream references(dataDirectory() + "/fitting/scipy_reference.csv");
    std::getline(references, line);
    while (std::getline(references, line)) {
        auto fields = splitCsvLine(line);
        auto &fitCase = cases[fields.at(0)];
        fitCase.distribution = fields.at(1);
        fitCase.parameters.emplace(fields.at(2), parseNumber(fields.at(3)));
    }
    return cases;
}

} // namespace

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
    const std::vector<double> nonFinite = {1, std::numeric_limits<double>::quiet_NaN()};
    CHECK_THROWS_AS(fit<NormalRand<double>>(nonFinite), fit_error);

    /// conjugate updates go through fit_bayes; fit must not return a likelihood estimate labelled Bayes
    CHECK_THROWS_WITH(fit<NormalRand<double>>(sample, Method::Bayes),
                      "Normal distribution: requested fitting method is not implemented");
    CHECK_THROWS_AS(fit<PoissonRand<int>>(std::vector<int>{1, 2, 3}, Method::Bayes), fit_error);
    CHECK(fit<PoissonRand<int>>(std::vector<int>{1, 2, 3}, Method::UMVU).distribution.GetRate() == 2);
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

    std::vector<double> fullLevySample;
    const LevyRand<double> levyReference(1.0, 2.0);
    for (int i = 0; i < 200; ++i)
        fullLevySample.push_back(levyReference.Quantile((i + 0.5) / 200));
    auto fullLevy = fit<LevyRand<double>>(fullLevySample);
    REQUIRE(fullLevy.converged);
    CHECK_THAT(fullLevy.distribution.GetLocation(), WithinAbs(1.0, 0.05));
    CHECK_THAT(fullLevy.distribution.GetScale(), WithinRel(2.0, 0.05));

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

TEST_CASE("Unified fitting covers Bernoulli and geometric distributions", "[fit]")
{
    const std::vector<int> bernoulliSample = {0, 1, 1, 0, 1};
    auto bernoulli = fit<BernoulliRand>(bernoulliSample);
    CHECK(bernoulli.converged);
    CHECK_THAT(bernoulli.distribution.GetProbability(), WithinAbs(0.6, 1e-14));

    const std::vector<int> geometricSample = {0, 1, 0, 2, 3, 1, 0};
    auto geometric = fit<GeometricRand<int>>(geometricSample);
    CHECK(geometric.converged);
    CHECK_THAT(geometric.distribution.GetProbability(), WithinAbs(7.0 / 14.0, 1e-14));
}

TEST_CASE("Generic MLE covers additional Tier 1 continuous distributions", "[fit]")
{
    auto quantileSample = []<class Distribution>(const Distribution &distribution) {
        constexpr int sampleSize = 200;
        std::vector<double> sample;
        sample.reserve(sampleSize);
        for (int i = 0; i < sampleSize; ++i)
            sample.push_back(distribution.Quantile((i + 0.5) / sampleSize));
        return sample;
    };

    auto cauchy = fit<CauchyRand<double>>(quantileSample(CauchyRand<double>(1.0, 2.0)));
    REQUIRE(cauchy.converged);
    CHECK_THAT(cauchy.distribution.GetLocation(), WithinAbs(1.0, 5e-8));
    CHECK_THAT(cauchy.distribution.GetScale(), WithinRel(2.0, 0.02));

    auto gumbel = fit<GumbelRand<double>>(quantileSample(GumbelRand<double>(-0.5, 1.5)));
    REQUIRE(gumbel.converged);
    CHECK_THAT(gumbel.distribution.GetLocation(), WithinAbs(-0.5, 0.02));
    CHECK_THAT(gumbel.distribution.GetScale(), WithinRel(1.5, 0.02));

    auto inverseGamma = fit<InverseGammaRand<double>>(quantileSample(InverseGammaRand<double>(4.0, 3.0)));
    REQUIRE(inverseGamma.converged);
    CHECK_THAT(inverseGamma.distribution.GetShape(), WithinRel(4.0, 0.03));
    CHECK_THAT(inverseGamma.distribution.GetRate(), WithinRel(3.0, 0.03));

    auto inverseGaussian = fit<InverseGaussianRand<double>>(quantileSample(InverseGaussianRand<double>(2.0, 5.0)));
    REQUIRE(inverseGaussian.converged);
    CHECK_THAT(inverseGaussian.distribution.GetMean(), WithinRel(2.0, 0.02));
    CHECK_THAT(inverseGaussian.distribution.GetShape(), WithinRel(5.0, 0.03));
}

TEST_CASE("Generic MLE supports free shape and noncentrality parameters", "[fit]")
{
    constexpr int sampleSize = 300;
    std::vector<double> studentSample;
    std::vector<double> chiSquaredSample;
    studentSample.reserve(sampleSize);
    chiSquaredSample.reserve(sampleSize);
    const StudentTRand<double> studentReference(7.0, 0.5, 1.2);
    const NoncentralChiSquaredRand<double> chiSquaredReference(4.0, 2.0);
    for (int i = 0; i < sampleSize; ++i) {
        const double probability = (i + 0.5) / sampleSize;
        studentSample.push_back(studentReference.Quantile(probability));
        chiSquaredSample.push_back(chiSquaredReference.Quantile(probability));
    }

    auto student = fit<StudentTRand<double>>(studentSample);
    REQUIRE(student.converged);
    CHECK_THAT(student.distribution.GetDegree(), WithinRel(7.0, 0.1));
    CHECK_THAT(student.distribution.GetLocation(), WithinAbs(0.5, 0.02));
    CHECK_THAT(student.distribution.GetScale(), WithinRel(1.2, 0.03));

    auto chiSquared = fit<NoncentralChiSquaredRand<double>>(chiSquaredSample);
    REQUIRE(chiSquared.converged);
    CHECK_THAT(chiSquared.distribution.GetDegree(), WithinRel(4.0, 0.05));
    CHECK_THAT(chiSquared.distribution.GetNoncentrality(), WithinRel(2.0, 0.1));
}

TEST_CASE("Generic MLE fits noncentral t noncentrality", "[fit][noncentral-t]")
{
    constexpr int sampleSize = 120;
    const NoncentralTRand<double> reference(8, 1.5);
    std::vector<double> sample;
    sample.reserve(sampleSize);
    for (int i = 0; i < sampleSize; ++i)
        sample.push_back(reference.Quantile((i + 0.5) / sampleSize));

    FitOptions options;
    options.fixed = {{"degree", 8}};
    auto result = fit<NoncentralTRand<double>>(sample, options);
    REQUIRE(result.converged);
    CHECK(result.distribution.GetDegree() == 8);
    CHECK_THAT(result.distribution.GetNoncentrality(), WithinAbs(1.5, 0.02));
    CHECK(result.parameter_names == std::vector<std::string>{"degree", "noncentrality"});
    CHECK(result.fixed == std::vector<bool>{true, false});
}

TEST_CASE("Generic MLE fits Frechet location, scale and shape", "[fit]")
{
    constexpr int sampleSize = 300;
    std::vector<double> sample;
    sample.reserve(sampleSize);
    const FrechetRand<double> reference(3.0, 2.0, -1.0);
    for (int i = 0; i < sampleSize; ++i)
        sample.push_back(reference.Quantile((i + 0.5) / sampleSize));

    auto result = fit<FrechetRand<double>>(sample);
    REQUIRE(result.converged);
    CHECK_THAT(result.distribution.GetShape(), WithinRel(3.0, 0.05));
    CHECK_THAT(result.distribution.GetScale(), WithinRel(2.0, 0.05));
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(-1.0, 0.05));
}

TEST_CASE("Generic MLE fits generalized Pareto and extreme-value distributions", "[fit]")
{
    constexpr int sampleSize = 400;
    std::vector<double> paretoSample;
    std::vector<double> extremeValueSample;
    paretoSample.reserve(sampleSize);
    extremeValueSample.reserve(sampleSize);
    const GeneralizedParetoRand<double> paretoReference(0.2, 2.0);
    const GeneralizedExtremeValueRand<double> extremeValueReference(0.2, 1.5, -0.5);
    for (int i = 0; i < sampleSize; ++i) {
        const double probability = (i + 0.5) / sampleSize;
        paretoSample.push_back(paretoReference.Quantile(probability));
        extremeValueSample.push_back(extremeValueReference.Quantile(probability));
    }

    auto pareto = fit<GeneralizedParetoRand<double>>(paretoSample);
    REQUIRE(pareto.converged);
    CHECK_THAT(pareto.distribution.GetShape(), WithinAbs(0.2, 0.02));
    CHECK_THAT(pareto.distribution.GetScale(), WithinRel(2.0, 0.02));
    CHECK(pareto.distribution.GetLocation() == 0);

    auto extremeValue = fit<GeneralizedExtremeValueRand<double>>(extremeValueSample);
    REQUIRE(extremeValue.converged);
    CHECK_THAT(extremeValue.distribution.GetShape(), WithinAbs(0.2, 0.02));
    CHECK_THAT(extremeValue.distribution.GetScale(), WithinRel(1.5, 0.02));
    CHECK_THAT(extremeValue.distribution.GetLocation(), WithinAbs(-0.5, 0.02));
}

TEST_CASE("Tail-index estimators expose order statistics for diagnostics", "[fit][tail]")
{
    const std::vector<double> hillSample = {1, 2, 4, 8, 16};
    CHECK_THAT(hill_estimator(hillSample, 2), WithinAbs(1.5 * std::log(2.0), 1e-14));

    const auto plot = hill_plot(hillSample, 3);
    REQUIRE(plot.size() == 3);
    CHECK(plot[0].order == 1);
    CHECK(plot[0].threshold == 8);
    CHECK_THAT(plot[0].estimate, WithinAbs(std::log(2.0), 1e-14));
    CHECK(plot[2].order == 3);
    CHECK(plot[2].threshold == 2);

    const std::vector<double> pickandsSample = {1, 2, 3, 4, 5, 6, 7, 8};
    CHECK_THAT(pickands_estimator(pickandsSample, 2), WithinAbs(-1.0, 1e-14));

    CHECK_THROWS_AS(hill_estimator(hillSample, 0), std::invalid_argument);
    CHECK_THROWS_AS(hill_estimator(std::vector<double>{0, 1}, 1), std::invalid_argument);
    CHECK_THROWS_AS(pickands_estimator(pickandsSample, 3), std::invalid_argument);
    CHECK_THROWS_AS(pickands_estimator(std::vector<double>{1, 1, 1, 1}, 1), std::invalid_argument);
}

TEST_CASE("Peaks-over-threshold fitting returns a shifted generalized Pareto law", "[fit][tail]")
{
    constexpr int sampleSize = 400;
    constexpr double threshold = 5;
    const GeneralizedParetoRand<double> reference(0.2, 2.0, threshold);
    std::vector<double> sample = {1, 3, threshold};
    sample.reserve(sampleSize + sample.size());
    for (int i = 0; i < sampleSize; ++i)
        sample.push_back(reference.Quantile((i + 0.5) / sampleSize));

    const auto result = fit_peaks_over_threshold(sample, threshold);
    REQUIRE(result.converged);
    CHECK_THAT(result.distribution.GetShape(), WithinAbs(0.2, 0.02));
    CHECK_THAT(result.distribution.GetScale(), WithinRel(2.0, 0.02));
    CHECK(result.distribution.GetLocation() == threshold);
    CHECK(result.parameter_names == std::vector<std::string>{"shape", "scale"});
}

TEST_CASE("Censored likelihood combines density and probability contributions", "[fit][censored]")
{
    using Observation = CensoredObservation<double>;
    const std::vector<Observation> observations = {Observation::exact(-0.2), Observation::left(-0.5),
                                                   Observation::right(0.5), Observation::interval(-0.1, 0.7)};
    FitOptions options;
    options.fixed = {{"location", 0}, {"variance", 1}};
    const auto result = fit_censored<NormalRand<double>>(observations, options);
    const NormalRand<double> normal;
    const double expected = normal.logf(-0.2) + std::log(normal.F(-0.5)) + std::log(normal.S(0.5)) +
                            std::log(normal.F(0.7) - normal.F(-0.1));

    REQUIRE(result.converged);
    CHECK_THAT(result.log_likelihood, WithinAbs(expected, 1e-13));
    CHECK(result.distribution.GetLocation() == 0);
    CHECK_THAT(static_cast<double>(result.distribution.Variance()), WithinAbs(1.0, 1e-14));
}

TEST_CASE("Censored and truncated normal MLE recover their source law", "[fit][censored]")
{
    using Observation = CensoredObservation<double>;
    constexpr int sampleSize = 500;
    constexpr double censoringPoint = 2;
    constexpr double lower = -1;
    constexpr double upper = 2.5;
    const NormalRand<double> reference(0.5, 2.25);
    std::vector<Observation> censored;
    std::vector<double> truncated;
    censored.reserve(sampleSize);
    truncated.reserve(sampleSize);
    const double lowerProbability = reference.F(lower);
    const double probabilityWidth = reference.F(upper) - lowerProbability;
    for (int i = 0; i < sampleSize; ++i) {
        const double probability = (i + 0.5) / sampleSize;
        const double value = reference.Quantile(probability);
        censored.push_back(value > censoringPoint ? Observation::right(censoringPoint) : Observation::exact(value));
        truncated.push_back(reference.Quantile(lowerProbability + probability * probabilityWidth));
    }

    const auto censoredResult = fit_censored<NormalRand<double>>(censored);
    REQUIRE(censoredResult.converged);
    CHECK_THAT(censoredResult.distribution.GetLocation(), WithinAbs(0.5, 0.03));
    CHECK_THAT(static_cast<double>(censoredResult.distribution.Variance()), WithinRel(2.25, 0.04));

    const auto truncatedResult = fit_truncated<NormalRand<double>>(truncated, lower, upper);
    REQUIRE(truncatedResult.converged);
    CHECK_THAT(truncatedResult.distribution.GetLocation(), WithinAbs(0.5, 0.03));
    CHECK_THAT(static_cast<double>(truncatedResult.distribution.Variance()), WithinRel(2.25, 0.04));
}

TEST_CASE("Censored and truncated fitting validate observation bounds", "[fit][censored]")
{
    using Observation = CensoredObservation<double>;
    CHECK_THROWS_AS(fit_censored<NormalRand<double>>(std::vector<Observation>{}), fit_error);
    CHECK_THROWS_AS(fit_censored<NormalRand<double>>(std::vector<Observation>{Observation::interval(2, 1)}),
                    std::invalid_argument);
    CHECK_THROWS_AS(fit_truncated<NormalRand<double>>(std::vector<double>{0, 2}, -1, 1), fit_error);
    CHECK_THROWS_AS(fit_truncated<NormalRand<double>>(std::vector<double>{0}, 1, 1), std::invalid_argument);
}

TEST_CASE("Unified fitting searches integer-valued Tier 1 parameters", "[fit]")
{
    constexpr int sampleSize = 200;
    std::vector<double> chiSquaredSample;
    std::vector<double> erlangSample;
    std::vector<double> fisherSample;
    chiSquaredSample.reserve(sampleSize);
    erlangSample.reserve(sampleSize);
    fisherSample.reserve(sampleSize);
    const ChiSquaredRand<double> chiSquaredReference(7);
    const ErlangRand<double> erlangReference(4, 2.0);
    const FisherFRand<double> fisherReference(5, 9);
    for (int i = 0; i < sampleSize; ++i) {
        const double probability = (i + 0.5) / sampleSize;
        chiSquaredSample.push_back(chiSquaredReference.Quantile(probability));
        erlangSample.push_back(erlangReference.Quantile(probability));
        fisherSample.push_back(fisherReference.Quantile(probability));
    }

    auto chiSquared = fit<ChiSquaredRand<double>>(chiSquaredSample);
    REQUIRE(chiSquared.converged);
    CHECK(chiSquared.distribution.GetDegree() == 7);
    CHECK(std::isnan(chiSquared.std_errors[0]));

    auto erlang = fit<ErlangRand<double>>(erlangSample);
    REQUIRE(erlang.converged);
    CHECK(erlang.distribution.GetShape() == 4);
    CHECK_THAT(erlang.distribution.GetRate(), WithinRel(2.0, 0.02));
    CHECK(std::isnan(erlang.std_errors[0]));
    CHECK(std::isfinite(erlang.std_errors[1]));

    auto fisher = fit<FisherFRand<double>>(fisherSample);
    REQUIRE(fisher.converged);
    CHECK(fisher.distribution.GetFirstDegree() == 5);
    CHECK(fisher.distribution.GetSecondDegree() == 9);

    FitOptions options;
    options.fixed = {{"shape", 4}};
    auto fixedErlang = fit<ErlangRand<double>>(erlangSample, options);
    CHECK(fixedErlang.converged);
    CHECK(fixedErlang.distribution.GetShape() == 4);
}

TEST_CASE("Conjugate Bayesian fitting returns posterior distributions", "[fit][bayes]")
{
    const std::vector<double> normalSample = {1, 2, 3};
    const NormalInverseGammaRand<double> normalPrior(0, 1, 2, 3);
    auto normalPosterior = fit_bayes(normalSample, normalPrior);
    CHECK_THAT(normalPosterior.GetLocation(), WithinAbs(1.5, 1e-14));
    CHECK_THAT(normalPosterior.GetPrecision(), WithinAbs(4.0, 1e-14));
    CHECK_THAT(normalPosterior.GetShape(), WithinAbs(3.5, 1e-14));
    CHECK_THAT(normalPosterior.GetRate(), WithinAbs(5.5, 1e-14));

    auto normalLocationPosterior = fit_bayes<NormalRand<double>>(normalSample, NormalRand<double>(0, 1), 2.0);
    CHECK_THAT(normalLocationPosterior.GetLocation(), WithinAbs(1.2, 1e-14));
    CHECK_THAT(static_cast<double>(normalLocationPosterior.Variance()), WithinAbs(0.4, 1e-14));

    auto normalVariancePosterior = fit_bayes<NormalRand<double>>(normalSample, InverseGammaRand<double>(2, 3), 2.0);
    CHECK_THAT(normalVariancePosterior.GetShape(), WithinAbs(3.5, 1e-14));
    CHECK_THAT(normalVariancePosterior.GetRate(), WithinAbs(4.0, 1e-14));

    const std::vector<double> logNormalSample = {std::exp(1.0), std::exp(2.0), std::exp(3.0)};
    auto logNormalPosterior = fit_bayes<LogNormalRand<double>>(logNormalSample, normalPrior);
    CHECK_THAT(logNormalPosterior.GetLocation(), WithinAbs(normalPosterior.GetLocation(), 1e-14));
    CHECK_THAT(logNormalPosterior.GetPrecision(), WithinAbs(normalPosterior.GetPrecision(), 1e-14));
    CHECK_THAT(logNormalPosterior.GetShape(), WithinAbs(normalPosterior.GetShape(), 1e-14));
    CHECK_THAT(logNormalPosterior.GetRate(), WithinAbs(normalPosterior.GetRate(), 1e-14));

    auto logNormalLocationPosterior = fit_bayes<LogNormalRand<double>>(logNormalSample, NormalRand<double>(0, 1), 2.0);
    CHECK_THAT(logNormalLocationPosterior.GetLocation(), WithinAbs(1.2, 1e-14));
    CHECK_THAT(static_cast<double>(logNormalLocationPosterior.Variance()), WithinAbs(0.4, 1e-14));

    auto logNormalVariancePosterior =
        fit_bayes<LogNormalRand<double>>(logNormalSample, InverseGammaRand<double>(2, 3), 2.0);
    CHECK_THAT(logNormalVariancePosterior.GetShape(), WithinAbs(3.5, 1e-14));
    CHECK_THAT(logNormalVariancePosterior.GetRate(), WithinAbs(4.0, 1e-14));

    const GammaRand<double> gammaPrior(2, 3);
    auto exponentialPosterior = fit_bayes(std::vector<double>{0.5, 1.0, 1.5}, gammaPrior);
    CHECK_THAT(exponentialPosterior.GetShape(), WithinAbs(5.0, 1e-14));
    CHECK_THAT(exponentialPosterior.GetRate(), WithinAbs(6.0, 1e-14));

    auto poissonPosterior = fit_bayes(std::vector<int>{1, 2, 3}, gammaPrior);
    CHECK_THAT(poissonPosterior.GetShape(), WithinAbs(8.0, 1e-14));
    CHECK_THAT(poissonPosterior.GetRate(), WithinAbs(6.0, 1e-14));

    const BetaRand<double> betaPrior(1, 1);
    auto binomialPosterior = fit_bayes<BinomialRand<int>>(std::vector<int>{1, 2, 3}, betaPrior, 5);
    CHECK_THAT(binomialPosterior.GetAlpha(), WithinAbs(7.0, 1e-14));
    CHECK_THAT(binomialPosterior.GetBeta(), WithinAbs(10.0, 1e-14));

    auto negativeBinomialPosterior =
        fit_bayes<NegativeBinomialRand<int, double>>(std::vector<int>{1, 2, 3}, betaPrior, 2.0);
    CHECK_THAT(negativeBinomialPosterior.GetAlpha(), WithinAbs(7.0, 1e-14));
    CHECK_THAT(negativeBinomialPosterior.GetBeta(), WithinAbs(7.0, 1e-14));

    const InverseGammaRand<double> inverseGammaPrior(2, 3);
    auto weibullPosterior = fit_bayes<WeibullRand<double>>(std::vector<double>{1, 2}, inverseGammaPrior, 2.0);
    CHECK_THAT(weibullPosterior.GetShape(), WithinAbs(4.0, 1e-14));
    CHECK_THAT(weibullPosterior.GetRate(), WithinAbs(8.0, 1e-14));

    auto paretoPosterior = fit_bayes<ParetoRand<double>>(std::vector<double>{1, 2, 4}, GammaRand<double>(2, 3), 1.0);
    CHECK_THAT(paretoPosterior.GetShape(), WithinAbs(5.0, 1e-14));
    CHECK_THAT(paretoPosterior.GetRate(), WithinAbs(3.0 + std::log(8.0), 1e-14));

    auto uniformPosterior =
        fit_bayes<UniformRand<double>>(std::vector<double>{0.2, 0.6, 0.8}, ParetoRand<double>(2, 0.5), 0.0);
    CHECK_THAT(uniformPosterior.GetShape(), WithinAbs(5.0, 1e-14));
    CHECK_THAT(uniformPosterior.GetScale(), WithinAbs(0.8, 1e-14));

    CHECK_THROWS_AS(fit_bayes<NormalRand<double>>(std::vector<double>{}, normalPrior), fit_error);
}

TEST_CASE("Method-of-moments fitting has a named entry point", "[fit]")
{
    const std::vector<double> sample = {1, 2, 3, 4};
    auto result = fit_moments<NormalRand<double>>(sample);
    CHECK(result.method == Method::MethodOfMoments);
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(2.5, 1e-14));
    CHECK_THAT(static_cast<double>(result.distribution.Variance()), WithinAbs(1.25, 1e-14));
}

TEST_CASE("Stable fitting uses a McCulloch quantile estimate in S0", "[fit][stable]")
{
    using Stable = StableRand<double>;
    const Stable reference(1.4, 0.4, 2.0, -1.0, Stable::Parameterization::S0);
    Pcg64 engine(2027);
    std::vector<double> sample(5000);
    reference.Sample(engine, sample);

    auto result = fit<Stable>(sample, Method::Quantile);
    CHECK(result.converged);
    CHECK(result.method == Method::Quantile);
    CHECK(result.distribution.GetParameterization() == Stable::Parameterization::S0);
    CHECK(result.parameter_names == std::vector<std::string>{"exponent", "skewness", "scale", "location"});
    CHECK_THAT(result.distribution.GetExponent(), WithinAbs(1.4, 0.08));
    CHECK_THAT(result.distribution.GetSkewness(), WithinAbs(0.4, 0.15));
    CHECK_THAT(result.distribution.GetScale(), WithinAbs(2.0, 0.15));
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(-1.0, 0.2));

    const std::array gaussianSample{-2.0, -1.0, 0.0, 1.0, 2.0};
    const auto gaussian = fit<Stable>(gaussianSample, Method::Quantile);
    CHECK(gaussian.distribution.GetExponent() == 2.0);
    CHECK(gaussian.distribution.GetSkewness() == 0.0);
}

TEST_CASE("Stable MLE refines the McCulloch estimate", "[fit][stable]")
{
    using Stable = StableRand<double>;
    const Stable reference(1.5, 0.3, 1.2, 0.5, Stable::Parameterization::S0);
    Pcg64 engine(2028);
    std::vector<double> sample(120);
    reference.Sample(engine, sample);

    FitOptions options;
    options.fixed = {{"exponent", 1.5}, {"skewness", 0.3}};
    options.max_iterations = 500;
    auto result = fit<Stable>(sample, options);
    CHECK(result.converged);
    CHECK(result.method == Method::MLE);
    CHECK(result.fixed == std::vector<bool>{true, true, false, false});
    CHECK_THAT(result.distribution.GetScale(), WithinAbs(1.2, 0.25));
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(0.5, 0.25));
}

TEST_CASE("Stable Count MLE supports free and fixed parameters", "[fit][stable-count]")
{
    using Distribution = StableCountRand<double>;
    constexpr int sampleSize = 80;
    const Distribution reference(0.5, 1.75, -0.4);
    std::vector<double> sample;
    sample.reserve(sampleSize);
    for (int i = 0; i < sampleSize; ++i)
        sample.push_back(reference.Quantile((i + 0.5) / sampleSize));

    FitOptions fixedOptions;
    fixedOptions.fixed = {{"exponent", 0.5}, {"location", -0.4}};
    auto fixed = fit<Distribution>(sample, fixedOptions);
    REQUIRE(fixed.converged);
    CHECK(fixed.fixed == std::vector<bool>{true, false, true});
    CHECK_THAT(fixed.distribution.GetScale(), WithinRel(1.75, 0.02));

    FitOptions freeOptions;
    freeOptions.max_iterations = 300;
    auto free = fit<Distribution>(sample, freeOptions);
    REQUIRE(free.converged);
    CHECK_THAT(free.distribution.GetExponent(), WithinAbs(0.5, 0.08));
    CHECK_THAT(free.distribution.GetScale(), WithinRel(1.75, 0.15));
    CHECK_THAT(free.distribution.GetLocation(), WithinAbs(-0.4, 0.2));

    CHECK_THROWS_AS(fit<Distribution>(std::array{1.0, 1.0, 1.0}), fit_error);
}

TEST_CASE("Log-stable fitting applies the stable estimator to logarithms", "[fit][stable][log-stable]")
{
    using Distribution = LogStableRand<double>;
    const Distribution reference(1.4, -0.4, 0.75, 0.25, Distribution::Parameterization::S0);
    Pcg64 engine(2029);
    std::vector<double> sample(5000);
    reference.Sample(engine, sample);

    auto result = fit<Distribution>(sample, Method::Quantile);
    CHECK(result.converged);
    CHECK(result.distribution.GetParameterization() == Distribution::Parameterization::S0);
    CHECK_THAT(result.distribution.GetExponent(), WithinAbs(1.4, 0.08));
    CHECK_THAT(result.distribution.GetSkewness(), WithinAbs(-0.4, 0.15));
    CHECK_THAT(result.distribution.GetScale(), WithinAbs(0.75, 0.08));
    CHECK_THAT(result.distribution.GetLocation(), WithinAbs(0.25, 0.1));
    CHECK_THROWS_AS(fit<Distribution>(std::array{0.0, 1.0}, Method::Quantile), fit_error);
}

TEST_CASE("Model selection chooses the smallest information criterion", "[fit]")
{
    constexpr int sampleSize = 400;
    std::vector<double> sample;
    sample.reserve(sampleSize);
    const ExponentialRand<double> reference(2.0);
    for (int i = 0; i < sampleSize; ++i)
        sample.push_back(reference.Quantile((i + 0.5) / sampleSize));

    auto result = fit_best<ExponentialRand<double>, GammaRand<double>, LogNormalRand<double>>(sample, Criterion::AIC);
    CHECK(result.index == 0);
    CHECK(result.criterion == Criterion::AIC);
    CHECK(std::holds_alternative<FitResult<ExponentialRand<double>>>(result.fit));
    CHECK_THAT(result.value, WithinAbs(std::get<0>(result.fit).aic, 0));
}

TEST_CASE("Goodness-of-fit reports continuous EDF statistics", "[fit][gof]")
{
    constexpr int sampleSize = 100;
    std::vector<double> sample;
    sample.reserve(sampleSize);
    const UniformRand<double> distribution;
    for (int i = 0; i < sampleSize; ++i)
        sample.push_back((i + 0.5) / sampleSize);

    auto result = goodness_of_fit(distribution, sample);
    CHECK_THAT(result.kolmogorov_smirnov, WithinAbs(0.005, 1e-14));
    CHECK(result.kolmogorov_smirnov_p_value > 0.99);
    CHECK(result.anderson_darling >= 0);
    CHECK(result.anderson_darling_p_value > 0.99);
    CHECK_THAT(result.cramer_von_mises, WithinAbs(1.0 / 1200.0, 1e-14));
    CHECK(result.cramer_von_mises_p_value > 0.99);

    const std::vector<double> referenceSample = {0.05, 0.1, 0.2, 0.35, 0.8, 0.9};
    auto reference = goodness_of_fit(distribution, referenceSample);
    CHECK_THAT(reference.cramer_von_mises, WithinAbs(0.1416666666666667, 1e-14));
    CHECK_THAT(reference.cramer_von_mises_p_value, WithinAbs(0.42549025501331084, 1e-10));
}

TEST_CASE("Goodness-of-fit reports a discrete chi-squared test", "[fit][gof]")
{
    std::vector<int> sample(100);
    std::fill(sample.begin() + 50, sample.end(), 1);
    auto result = goodness_of_fit(BernoulliRand(0.5), sample);

    CHECK_THAT(result.chi_squared, WithinAbs(0, 1e-14));
    CHECK_THAT(result.p_value, WithinAbs(1, 1e-14));
    CHECK(result.degrees_of_freedom == 1);
    CHECK(result.bins == 2);

    CHECK_THROWS_AS(goodness_of_fit(BernoulliRand(0.5), std::vector<int>{0, 1}), std::invalid_argument);
}

TEST_CASE("Member goodness-of-fit tests validate samples and count every cell", "[fit][gof]")
{
    /// Binomial(4, 1/2) expects 10, 40, 60, 40 and 10 of 160 observations
    const BinomialRand<int> binomial(4, 0.5);
    auto sample = [](std::initializer_list<int> counts) {
        std::vector<int> values;
        int value = 0;
        for (int count : counts) {
            values.insert(values.end(), count, value);
            ++value;
        }
        return values;
    };
    /// Pearson's statistic 8 lies between the 95% quantiles of chi-squared with 3 and 4 degrees of freedom
    const auto shifted = sample({18, 32, 60, 40, 10});
    CHECK(binomial.PearsonChiSquaredTest(shifted, 0.05));
    CHECK_FALSE(binomial.PearsonChiSquaredTest(shifted, 0.05, 1));
    /// the empty cell {4} contributes its expected count 10 to the statistic 12.5
    CHECK_FALSE(binomial.PearsonChiSquaredTest(sample({10, 40, 60, 50}), 0.05));

    const std::vector<int> unsorted = {2, 1, 2, 3};
    CHECK_THROWS_AS(binomial.PearsonChiSquaredTest(unsorted, 0.05), std::invalid_argument);
    CHECK_THROWS_AS(binomial.PearsonChiSquaredTest(std::span<const int>{}, 0.05), std::invalid_argument);
    CHECK_THROWS_AS(binomial.PearsonChiSquaredTest(shifted, 1.5), std::invalid_argument);
    /// every observation in the lower cell (-∞, 2]: the statistic 5.45 has four degrees of freedom
    const PoissonRand<int> poisson(3);
    bool accepted = false;
    REQUIRE_NOTHROW(accepted = poisson.PearsonChiSquaredTest(std::vector<int>{0, 0, 1, 1}, 0.05, 2, 5));
    CHECK(accepted);
    /// unobserved values of a heavy tail are accounted for together, not one by one
    ZetaRand<int> zeta(1.2);
    zeta.Reseed(5);
    std::vector<int> heavy(2000);
    zeta.Sample(heavy);
    std::ranges::sort(heavy);
    REQUIRE(heavy.back() > 1000000);
    CHECK_NOTHROW(zeta.PearsonChiSquaredTest(heavy, 0.05));

    const NormalRand<double> normal;
    CHECK_THROWS_AS(normal.KolmogorovSmirnovTest(std::span<const double>{}, 0.05), std::invalid_argument);
    CHECK_THROWS_AS(normal.KolmogorovSmirnovTest(std::vector<double>{0.0}, 0.0), std::invalid_argument);
}

TEST_CASE("Unified fitting agrees with SciPy on fixed datasets", "[fit][scipy]")
{
    const auto cases = loadScipyFitCases();
    auto check = [&]<class Distribution>(std::string_view name) {
        const auto &fitCase = cases.at(std::string(name));
        INFO(name << " / " << fitCase.distribution);
        auto result = fit<Distribution>(fitCase.sample);
        REQUIRE(result.converged);
        REQUIRE(result.parameter_names.size() == result.parameters.size());
        for (std::size_t i = 0; i < result.parameters.size(); ++i) {
            const double expected = fitCase.parameters.at(result.parameter_names[i]);
            CHECK_THAT(result.parameters[i], WithinAbs(expected, 5e-5 * std::max(1.0, std::fabs(expected))));
        }
    };

    check.template operator()<NormalRand<double>>("normal");
    check.template operator()<LogNormalRand<double>>("log_normal");
    check.template operator()<GammaRand<double>>("gamma");
    check.template operator()<BetaRand<double>>("beta");
    check.template operator()<WeibullRand<double>>("weibull");
    check.template operator()<CauchyRand<double>>("cauchy");
    check.template operator()<GumbelRand<double>>("gumbel");
    check.template operator()<InverseGammaRand<double>>("inverse_gamma");
    check.template operator()<LogisticRand<double>>("logistic");
}

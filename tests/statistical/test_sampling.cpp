// Goodness of fit and moments of the samplers, with fixed seeds so that results are
// deterministic. If one of these fails after a sampler change, investigate the sampler:
// never change the seed until the test passes.

#include <catch2/catch_test_macros.hpp>

#include "registry.h"
#include "stats.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <map>
#include <numeric>
#include <sstream>

using namespace randlib_test;

namespace
{

constexpr size_t SAMPLE_SIZE = 200000;
constexpr size_t MIN_SAMPLE_SIZE = 20000;
constexpr double TIME_BUDGET_SECONDS = 4.0;
constexpr double MIN_P_VALUE = 1e-4;
constexpr double MAX_STANDARD_ERRORS = 5.0;
constexpr unsigned long BASE_SEED = 20231010;

/// Sample size within the time budget for evaluating the cdf once per draw.
template <class Cdf>
size_t affordableSampleSize(const Cdf &cdf, const std::vector<double> &probe)
{
    auto start = std::chrono::steady_clock::now();
    for (double x : probe)
        cdf(x);
    double perCall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() / probe.size();
    size_t affordable = static_cast<size_t>(TIME_BUDGET_SECONDS / std::max(perCall, 1e-9));
    return std::clamp(affordable, MIN_SAMPLE_SIZE, SAMPLE_SIZE);
}

template <class T>
std::string checkMoments(const UnivariateDistribution<T> &d, const std::vector<T> &sample)
{
    std::ostringstream out;
    double n = static_cast<double>(sample.size());
    double mean = static_cast<double>(d.Mean()), variance = static_cast<double>(d.Variance());
    if (!std::isfinite(mean) || !std::isfinite(variance))
        return {};
    long double sum = 0, sumSq = 0;
    for (T x : sample)
        sum += x;
    long double sampleMean = sum / n;
    for (T x : sample)
        sumSq += (x - sampleMean) * (x - sampleMean);
    double sampleVariance = static_cast<double>(sumSq / (n - 1));

    double meanError = std::fabs(static_cast<double>(sampleMean) - mean) / std::sqrt(variance / n);
    if (meanError > MAX_STANDARD_ERRORS)
        out << "sample mean " << static_cast<double>(sampleMean) << " is " << meanError
            << " standard errors from Mean() = " << mean << "; ";

    double kurtosis = static_cast<double>(d.ExcessKurtosis());
    if (std::isfinite(kurtosis)) {
        /// Var(s²) = σ⁴ (κ / n + 2 / (n - 1)) exactly, κ being the excess kurtosis
        double varianceError =
            std::fabs(sampleVariance - variance) / (variance * std::sqrt(kurtosis / n + 2 / (n - 1)));
        if (varianceError > MAX_STANDARD_ERRORS)
            out << "sample variance " << sampleVariance << " is " << varianceError
                << " standard errors from Variance() = " << variance << "; ";
    }
    return out.str();
}

/// Empirical characteristic function against CF(t), at frequencies scaled to the sample spread.
/// |mean(exp(itX)) - φ(t)| has standard deviation at most 1/√n.
template <class T>
std::string checkCharacteristicFunction(const UnivariateDistribution<T> &d, std::vector<double> sorted)
{
    std::ostringstream out;
    double n = static_cast<double>(sorted.size());
    double spread = sorted[static_cast<size_t>(0.75 * (n - 1))] - sorted[static_cast<size_t>(0.25 * (n - 1))];
    if (!(spread > 0))
        return {};
    for (double scaledT : {0.3, 1.0, 3.0}) {
        double t = scaledT / spread;
        std::complex<double> phi;
        try {
            phi = d.CF(t);
        } catch (const randlib::not_implemented &) {
            if (d.Name().starts_with("Log-Stable("))
                return {};
            return "characteristic function unexpectedly throws not_implemented; ";
        }
        long double re = 0, im = 0;
        for (double x : sorted) {
            re += std::cos(t * x);
            im += std::sin(t * x);
        }
        std::complex<double> empirical(static_cast<double>(re / n), static_cast<double>(im / n));
        double distance = std::abs(empirical - phi) * std::sqrt(n);
        if (distance > MAX_STANDARD_ERRORS)
            out << "empirical characteristic function at t = " << t << " is " << empirical << ", CF(t) = " << phi
                << "; ";
    }
    return out.str();
}

template <class Make>
std::vector<std::string> runContinuous(const std::vector<const Case *> &cases, Make make)
{
    std::vector<std::string> failures(cases.size());
    parallelFor(cases.size(), [&](size_t i) {
        std::ostringstream out;
        try {
            auto d = make(*cases[i]);
            d->Reseed(BASE_SEED + i);
            std::vector<double> probe(200);
            d->Sample(probe);
            bool hasCdf = true;
            try {
                static_cast<void>(d->F(probe.front()));
            } catch (const randlib::not_implemented &) {
                hasCdf = false;
            }
            size_t n = hasCdf ? affordableSampleSize([&d](double x) { return d->F(x); }, probe) : SAMPLE_SIZE;
            std::vector<double> sample(n);
            d->Sample(sample);
            if (std::any_of(sample.begin(), sample.end(),
                            [&d](double x) { return x < d->MinValue() || x > d->MaxValue() || std::isnan(x); }))
                out << "sample outside the support; ";
            std::string moments = checkMoments(*d, sample);
            std::sort(sample.begin(), sample.end());
            if (hasCdf) {
                double D = kolmogorovSmirnovStatistic(sample, [&d](double x) { return d->F(x); });
                double p = kolmogorovSmirnovPValue(D, n);
                if (p < MIN_P_VALUE)
                    out << "Kolmogorov-Smirnov p-value " << p << " (D = " << D << ", n = " << n << "); ";
            }
            out << moments << checkCharacteristicFunction(*d, sample);
        } catch (const std::exception &e) {
            out << "exception: " << e.what();
        }
        failures[i] = out.str();
    });
    return failures;
}

void report(const std::vector<const Case *> &cases, const std::vector<std::string> &failures)
{
    size_t total = 0;
    for (size_t i = 0; i < cases.size(); ++i) {
        if (!failures[i].empty()) {
            FAIL_CHECK(cases[i]->id + ": " + failures[i]);
            ++total;
        }
    }
    CHECK(total == 0);
}

std::string caseId(std::string type, const std::string &args)
{
    const std::string plain = "RL_PLAIN(";
    if (type.rfind(plain, 0) == 0)
        type = type.substr(plain.size(), type.size() - plain.size() - 1);
    auto lt = type.find('<');
    if (lt != std::string::npos)
        type.resize(lt);
    return type + "(" + args + ")";
}

void checkDiscreteCase(const std::string &type, const std::string &args)
{
    const std::string id = caseId(type, args);
    const auto cases = casesOfKind(Kind::Discrete);
    const auto selected = std::find_if(cases.begin(), cases.end(), [&id](const Case *c) { return c->id == id; });
    REQUIRE(selected != cases.end());
    const size_t index = static_cast<size_t>(std::distance(cases.begin(), selected));
    std::ostringstream out;
    try {
        auto d = (*selected)->makeDiscrete();
        d->Reseed(BASE_SEED + index);
        std::vector<int> sample(SAMPLE_SIZE);
        d->Sample(sample);
        if (std::any_of(sample.begin(), sample.end(), [&d](int k) { return k < d->MinValue() || k > d->MaxValue(); }))
            out << "sample outside the support; ";
        out << checkMoments(*d, sample);

        std::map<int, size_t> counts;
        for (int k : sample)
            ++counts[k];
        const double n = static_cast<double>(sample.size());
        std::vector<std::pair<double, double>> bins;
        double observed = 0, expected = n * d->F(counts.begin()->first - 1);
        for (auto it = counts.begin(); it != counts.end(); ++it) {
            observed += static_cast<double>(it->second);
            expected += n * d->P(it->first);
            auto next = std::next(it);
            int gapEnd = (next == counts.end()) ? it->first : next->first - 1;
            if (gapEnd > it->first)
                expected += n * (d->F(gapEnd) - d->F(it->first));
            if (next == counts.end())
                expected += n * d->S(it->first);
            if (expected >= 5) {
                bins.push_back({observed, expected});
                observed = expected = 0;
            }
        }
        if (!bins.empty()) {
            bins.back().first += observed;
            bins.back().second += expected;
        }
        if (bins.size() > 1) {
            double chiSq = 0;
            for (auto [o, e] : bins)
                chiSq += (o - e) * (o - e) / e;
            double p = chiSquaredSurvival(chiSq, static_cast<double>(bins.size() - 1));
            if (p < MIN_P_VALUE)
                out << "chi-squared p-value " << p << " (statistic " << chiSq << ", " << bins.size() << " bins); ";
        }
    } catch (const std::exception &e) {
        out << "exception: " << e.what();
    }
    INFO(id);
    CHECK(out.str().empty());
}

} // namespace

TEST_CASE("Continuous samples follow their distribution", "[statistical]")
{
    auto cases = casesOfKind(Kind::Continuous);
    std::erase_if(cases, [](const Case *testCase) { return !testCase->testSampling; });
    report(cases, runContinuous(cases, [](const Case &c) { return c.makeContinuous(); }));
}

TEST_CASE("Stable Count inverse samples match independent quartiles", "[statistical][stable-count]")
{
    struct SamplingCase
    {
        double exponent;
        std::array<double, 3> quartiles;
    };
    constexpr std::array cases{
        SamplingCase{0.2, {981.1353613230349, 3865.491718469079, 13278.848675811233}},
        SamplingCase{0.8, {1.0920923194471634, 1.4831845398887014, 1.8506681858753713}},
        SamplingCase{0.95, {1.00947313942679, 1.1232778646876066, 1.199141995676045}},
    };
    constexpr std::size_t sampleSize = 1024;
    for (std::size_t caseIndex = 0; caseIndex < cases.size(); ++caseIndex) {
        const auto &[exponent, quartiles] = cases[caseIndex];
        StableCountRand<double> distribution(exponent);
        Pcg64 engine(BASE_SEED + caseIndex);
        std::array<std::size_t, 4> counts{};
        for (std::size_t i = 0; i < sampleSize; ++i) {
            const double value = distribution(engine);
            const auto bin = std::upper_bound(quartiles.begin(), quartiles.end(), value) - quartiles.begin();
            ++counts[static_cast<std::size_t>(bin)];
        }
        double statistic = 0;
        for (std::size_t count : counts) {
            constexpr double expected = sampleSize / 4.0;
            statistic += (count - expected) * (count - expected) / expected;
        }
        INFO("alpha = " << exponent);
        CHECK(chiSquaredSurvival(statistic, 3) > MIN_P_VALUE);
    }
}

TEST_CASE("Static binomial variates follow the binomial law", "[statistical][discrete]")
{
    struct BinomialCase
    {
        int number;
        double probability;
    };
    constexpr std::array cases{BinomialCase{12, 0.3},  BinomialCase{50, 0.15},      BinomialCase{50, 0.85},
                               BinomialCase{200, 0.3}, BinomialCase{1000000, 1e-5}, BinomialCase{1000000, 0.4}};
    for (std::size_t caseIndex = 0; caseIndex < cases.size(); ++caseIndex) {
        const auto [number, probability] = cases[caseIndex];
        const BinomialRand<int> distribution(number, probability);
        Pcg64 engine(BASE_SEED + caseIndex);
        std::map<int, std::size_t> counts;
        for (std::size_t i = 0; i < SAMPLE_SIZE; ++i)
            ++counts[BinomialRand<int>::Variate(number, probability, engine)];
        double statistic = 0, expected = distribution.F(counts.begin()->first - 1) * SAMPLE_SIZE, observed = 0;
        int cells = 0;
        for (int k = counts.begin()->first; k <= counts.rbegin()->first; ++k) {
            expected += distribution.P(k) * SAMPLE_SIZE;
            observed += counts.contains(k) ? counts[k] : 0;
            if (expected >= 20) {
                statistic += (observed - expected) * (observed - expected) / expected;
                ++cells;
                expected = observed = 0;
            }
        }
        expected += distribution.S(counts.rbegin()->first) * SAMPLE_SIZE;
        statistic += (observed - expected) * (observed - expected) / std::max(expected, 1.0);
        INFO("Binomial(" << number << ", " << probability << ")");
        CHECK(chiSquaredSurvival(statistic, cells) > MIN_P_VALUE);
    }
}

TEST_CASE("Samplers with large parameters keep bounded work", "[statistical][discrete]")
{
    /// each of these took time linear in a parameter of order 10^8 per draw before
    constexpr std::size_t sampleSize = 4000;
    auto check = [](const std::string &name, const auto &distribution, unsigned long seed) {
        Pcg64 engine(seed);
        std::vector<long long> sample(sampleSize);
        for (auto &x : sample)
            x = distribution(engine);
        const double mean = std::accumulate(sample.begin(), sample.end(), 0.0) / sampleSize;
        const double standardError = std::sqrt(static_cast<double>(distribution.Variance()) / sampleSize);
        INFO(name << ": sample mean " << mean << ", Mean() " << static_cast<double>(distribution.Mean()));
        CHECK(std::fabs(mean - static_cast<double>(distribution.Mean())) < MAX_STANDARD_ERRORS * standardError);
    };
    check("Pascal(100000000, 0.5)", PascalRand<long long>(100000000, 0.5), BASE_SEED);
    check("HyperGeometric(1000000000, 300000000, 600000000)",
          HyperGeometricRand<long long>(1000000000, 300000000, 600000000), BASE_SEED + 1);
    check("NegativeHyperGeometric(1000000000, 500000000, 200000000)",
          NegativeHyperGeometricRand<long long>(1000000000, 500000000, 200000000), BASE_SEED + 2);
    check("BetaBinomial(1000000000, 3, 4)", BetaBinomialRand<long long>(1000000000, 3, 4), BASE_SEED + 3);
}

TEST_CASE("Singular samples follow their distribution", "[statistical]")
{
    auto cases = casesOfKind(Kind::Singular);
    report(cases, runContinuous(cases, [](const Case &c) { return c.makeSingular(); }));
}

#define RL_CONT(T, ...)
#define RL_CONT_V2(T, ...)
#define RL_CONT_V2_NO_STATS(T, ...)
#define RL_SING(T, ...)
#define RL_BIVC(T, ...)
#define RL_BIVD(T, ...)
#define RL_DISC(T, ...)                                                                                                \
    TEST_CASE("Discrete samples follow " #T "(" #__VA_ARGS__ ")", "[statistical][discrete]")                           \
    {                                                                                                                  \
        checkDiscreteCase(#T, #__VA_ARGS__);                                                                           \
    }

#include "../support/cases.inc"

#undef RL_CONT
#undef RL_CONT_V2
#undef RL_CONT_V2_NO_STATS
#undef RL_DISC
#undef RL_SING
#undef RL_BIVC
#undef RL_BIVD

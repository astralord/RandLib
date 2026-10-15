#pragma once

#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/kolmogorov_smirnov.hpp>
#include <randlib/math/special_functions.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <map>
#include <numbers>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace randlib
{

struct ContinuousGoodnessOfFit
{
    double kolmogorov_smirnov = 0;
    double kolmogorov_smirnov_p_value = 1;
    double anderson_darling = 0;
    double anderson_darling_p_value = 1;
    double cramer_von_mises = 0;
    double cramer_von_mises_p_value = 1;
};

struct DiscreteGoodnessOfFit
{
    double chi_squared = 0;
    double p_value = 1;
    int degrees_of_freedom = 0;
    std::size_t bins = 0;
};

namespace detail
{

inline double andersonDarlingCdf(double statistic)
{
    if (!(statistic > 0))
        return 0;
    if (statistic < 2) {
        const double polynomial =
            2.00012 +
            statistic *
                (0.247105 -
                 statistic * (0.0649821 - statistic * (0.0347962 - statistic * (0.011672 - 0.00168691 * statistic))));
        return std::clamp(std::exp(-1.2337141 / statistic) * polynomial / std::sqrt(statistic), 0.0, 1.0);
    }
    const double exponent =
        1.0776 -
        statistic *
            (2.30695 - statistic * (0.43424 - statistic * (0.082433 - statistic * (0.008056 - 0.0003146 * statistic))));
    return std::clamp(std::exp(-std::exp(exponent)), 0.0, 1.0);
}

inline double cvmCdfInfinite(double statistic)
{
    if (!(statistic > 0))
        return 0;
    if (statistic >= 4)
        return 1;
    double sum = 0;
    for (int k = 0; k < 100; ++k) {
        const double y = 4 * k + 1;
        const double q = y * y / (16 * statistic);
        const double gammaRatio = std::exp(std::lgamma(k + 0.5) - std::lgamma(k + 1.0));
        const double term = gammaRatio * std::sqrt(y) * std::exp(-q) * math::besselK(0.25, q) /
                            (std::pow(std::numbers::pi, 1.5) * std::sqrt(statistic));
        sum += term;
        if (std::fabs(term) < 1e-12)
            break;
    }
    return std::clamp(sum, 0.0, 1.0);
}

inline double cvmEd2(double value)
{
    const double z = value * value / 4;
    return std::exp(-z) * std::pow(value / 2, 1.5) * (math::besselK(0.25, z) + math::besselK(0.75, z)) /
           std::sqrt(std::numbers::pi);
}

inline double cvmEd3(double value)
{
    const double z = value * value / 4;
    const double bessel = 2 * math::besselK(0.25, z) + 3 * math::besselK(0.75, z) - math::besselK(1.25, z);
    return std::exp(-z) * std::pow(value / 2, 2.5) * bessel / std::sqrt(std::numbers::pi);
}

inline double cvmCorrectionTerm(int k, double statistic)
{
    const double m = 2 * k + 1;
    const double sx = 2 * std::sqrt(statistic);
    const double y1 = std::pow(statistic, 0.75);
    const double y2 = std::pow(statistic, 1.25);
    const double gamma1 = std::tgamma(k + 0.5);
    const double gamma2 = std::tgamma(k + 1.5);
    const double e1 = m * gamma1 * cvmEd2((4 * k + 3) / sx) / (9 * y1);
    const double e2 = gamma1 * cvmEd3((4 * k + 1) / sx) / (72 * y2);
    const double e3 = 2 * (m + 2) * gamma2 * cvmEd3((4 * k + 5) / sx) / (12 * y2);
    const double e4 = 7 * m * gamma1 * cvmEd2((4 * k + 1) / sx) / (144 * y1);
    const double e5 = 7 * m * gamma1 * cvmEd2((4 * k + 5) / sx) / (144 * y1);
    return e1 + e2 + e3 + e4 + e5;
}

inline double cvmCdf(double statistic, std::size_t sampleSize)
{
    const double minimum = 1.0 / (12 * sampleSize);
    const double maximum = static_cast<double>(sampleSize) / 3;
    if (statistic <= minimum)
        return 0;
    if (statistic >= maximum)
        return 1;
    const double infinite = cvmCdfInfinite(statistic);
    double correction = 0;
    for (int k = 0; k < 100; ++k) {
        const double term = -cvmCorrectionTerm(k, statistic) / (std::numbers::pi * std::tgamma(k + 1.0));
        correction += term;
        if (std::fabs(term) < 1e-10)
            break;
    }
    return std::clamp(infinite * (1 + 1.0 / (12 * sampleSize)) + correction / sampleSize, 0.0, 1.0);
}

struct ChiSquaredBin
{
    double observed;
    double expected;
};

inline void mergeSparseBins(std::vector<ChiSquaredBin> &bins)
{
    while (bins.size() > 2) {
        const auto sparse = std::ranges::min_element(bins, {}, &ChiSquaredBin::expected);
        if (sparse->expected >= 5)
            return;
        const std::size_t index = static_cast<std::size_t>(sparse - bins.begin());
        const std::size_t target = index == 0 ? 1 : index - 1;
        bins[target].observed += sparse->observed;
        bins[target].expected += sparse->expected;
        bins.erase(sparse);
    }
}

} // namespace detail

template <class Distribution>
    requires std::floating_point<typename Distribution::value_type>
ContinuousGoodnessOfFit goodness_of_fit(const Distribution &distribution,
                                        std::span<const typename Distribution::value_type> sample)
{
    if (sample.empty())
        throw std::invalid_argument("Goodness-of-fit testing requires a non-empty sample");
    std::vector<typename Distribution::value_type> sorted(sample.begin(), sample.end());
    std::ranges::sort(sorted);
    const double n = static_cast<double>(sorted.size());
    double kolmogorovSmirnov = 0;
    double andersonDarlingSum = 0;
    double cramerVonMises = 1 / (12 * n);
    const double minimumProbability = std::numeric_limits<double>::min();
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        const double rawProbability = distribution.F(sorted[i]);
        const double rawReverseProbability = distribution.F(sorted[sorted.size() - i - 1]);
        if (!std::isfinite(rawProbability) || !std::isfinite(rawReverseProbability) || rawProbability < 0 ||
            rawProbability > 1 || rawReverseProbability < 0 || rawReverseProbability > 1)
            throw std::runtime_error("The fitted distribution returned an invalid probability");
        const double probability = std::clamp(rawProbability, minimumProbability, 1.0);
        const double reverseProbability =
            std::clamp(rawReverseProbability, 0.0, 1 - std::numeric_limits<double>::epsilon());
        kolmogorovSmirnov = std::max(kolmogorovSmirnov, (i + 1) / n - probability);
        kolmogorovSmirnov = std::max(kolmogorovSmirnov, probability - i / n);
        andersonDarlingSum += (2 * i + 1) * (std::log(probability) + std::log1p(-reverseProbability));
        const double expected = (2 * i + 1) / (2 * n);
        cramerVonMises += (probability - expected) * (probability - expected);
    }

    ContinuousGoodnessOfFit result;
    result.kolmogorov_smirnov = kolmogorovSmirnov;
    result.kolmogorov_smirnov_p_value = KolmogorovSmirnovRand<double>().S(std::sqrt(n) * kolmogorovSmirnov);
    result.anderson_darling = -n - andersonDarlingSum / n;
    result.anderson_darling_p_value = 1 - detail::andersonDarlingCdf(result.anderson_darling);
    result.cramer_von_mises = cramerVonMises;
    result.cramer_von_mises_p_value = 1 - detail::cvmCdf(cramerVonMises, sorted.size());
    return result;
}

template <class Distribution>
    requires std::integral<typename Distribution::value_type>
DiscreteGoodnessOfFit goodness_of_fit(const Distribution &distribution,
                                      std::span<const typename Distribution::value_type> sample,
                                      int estimatedParameters = 0)
{
    if (sample.empty())
        throw std::invalid_argument("Goodness-of-fit testing requires a non-empty sample");
    if (estimatedParameters < 0)
        throw std::invalid_argument("The number of estimated parameters should be non-negative");
    std::map<typename Distribution::value_type, std::size_t> counts;
    for (auto value : sample)
        ++counts[value];
    std::vector<detail::ChiSquaredBin> bins;
    bins.reserve(counts.size() + 1);
    double totalExpected = 0;
    for (const auto &[value, observed] : counts) {
        const double probability = distribution.P(value);
        if (!std::isfinite(probability) || probability < 0 || probability > 1)
            throw std::runtime_error("The fitted distribution returned an invalid probability");
        const double expected = sample.size() * probability;
        bins.push_back({static_cast<double>(observed), expected});
        totalExpected += expected;
    }
    if (totalExpected > sample.size() * (1 + 16 * std::numeric_limits<double>::epsilon()))
        throw std::runtime_error("The fitted distribution returned probabilities summing above one");
    bins.push_back({0, std::max(0.0, sample.size() - totalExpected)});
    detail::mergeSparseBins(bins);
    if (std::ranges::any_of(bins, [](const auto &bin) { return bin.expected < 5; }))
        throw std::invalid_argument("Chi-squared testing requires expected counts of at least five");
    const int degreesOfFreedom = static_cast<int>(bins.size()) - 1 - estimatedParameters;
    if (degreesOfFreedom <= 0)
        throw std::invalid_argument("Too few chi-squared bins for the estimated parameter count");
    double statistic = 0;
    for (const auto &bin : bins) {
        if (!(bin.expected > 0))
            throw std::invalid_argument("The fitted distribution assigns zero probability to an observed bin");
        const double difference = bin.observed - bin.expected;
        statistic += difference * difference / bin.expected;
    }
    return {statistic, ChiSquaredRand<double>(degreesOfFreedom).S(statistic), degreesOfFreedom, bins.size()};
}

template <class Distribution, std::ranges::contiguous_range Range>
    requires std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, typename Distribution::value_type>
auto goodness_of_fit(const Distribution &distribution, const Range &sample)
{
    return goodness_of_fit(distribution, std::span<const typename Distribution::value_type>(sample));
}

template <class Distribution, std::ranges::contiguous_range Range>
    requires std::integral<typename Distribution::value_type> &&
             std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, typename Distribution::value_type>
DiscreteGoodnessOfFit goodness_of_fit(const Distribution &distribution, const Range &sample, int estimatedParameters)
{
    return goodness_of_fit(distribution, std::span<const typename Distribution::value_type>(sample),
                           estimatedParameters);
}

} // namespace randlib

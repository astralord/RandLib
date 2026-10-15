#pragma once

#include <randlib/continuous/generalized_pareto.hpp>
#include <randlib/fit/mle.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace randlib
{

struct TailIndexEstimate
{
    std::size_t order = 0;
    double threshold = 0;
    double estimate = 0;
};

namespace detail
{

template <std::floating_point Real>
std::vector<Real> sortedTailSample(std::span<const Real> sample)
{
    if (sample.size() < 2)
        throw std::invalid_argument("Tail-index estimation requires at least two observations");
    if (std::ranges::any_of(sample, [](Real value) { return !(value > 0) || !std::isfinite(value); }))
        throw std::invalid_argument("Tail-index estimation requires positive finite observations");
    std::vector<Real> sorted(sample.begin(), sample.end());
    std::ranges::sort(sorted);
    return sorted;
}

template <std::floating_point Real>
double hillFromSorted(std::span<const Real> sorted, std::size_t order)
{
    if (order == 0 || order >= sorted.size())
        throw std::invalid_argument("Hill estimator order should be between one and sample size minus one");
    const double logThreshold = std::log(static_cast<double>(sorted[sorted.size() - order - 1]));
    long double sum = 0;
    for (std::size_t i = sorted.size() - order; i < sorted.size(); ++i)
        sum += std::log(static_cast<double>(sorted[i])) - logThreshold;
    return static_cast<double>(sum / order);
}

} // namespace detail

template <std::floating_point Real>
double hill_estimator(std::span<const Real> sample, std::size_t order)
{
    const auto sorted = detail::sortedTailSample(sample);
    return detail::hillFromSorted(std::span<const Real>(sorted), order);
}

template <std::ranges::contiguous_range Range>
    requires std::floating_point<std::remove_cv_t<std::ranges::range_value_t<Range>>> &&
             (!std::same_as<std::remove_cvref_t<Range>,
                            std::span<const std::remove_cv_t<std::ranges::range_value_t<Range>>>>)
double hill_estimator(const Range &sample, std::size_t order)
{
    using Real = std::remove_cv_t<std::ranges::range_value_t<Range>>;
    return hill_estimator(std::span<const Real>(sample), order);
}

template <std::floating_point Real>
std::vector<TailIndexEstimate> hill_plot(std::span<const Real> sample, std::size_t maximumOrder = 0)
{
    const auto sorted = detail::sortedTailSample(sample);
    if (maximumOrder == 0)
        maximumOrder = sorted.size() - 1;
    if (maximumOrder >= sorted.size())
        throw std::invalid_argument("Maximum Hill order should be smaller than the sample size");

    std::vector<TailIndexEstimate> result;
    result.reserve(maximumOrder);
    for (std::size_t order = 1; order <= maximumOrder; ++order) {
        result.push_back({order, static_cast<double>(sorted[sorted.size() - order - 1]),
                          detail::hillFromSorted(std::span<const Real>(sorted), order)});
    }
    return result;
}

template <std::ranges::contiguous_range Range>
    requires std::floating_point<std::remove_cv_t<std::ranges::range_value_t<Range>>> &&
             (!std::same_as<std::remove_cvref_t<Range>,
                            std::span<const std::remove_cv_t<std::ranges::range_value_t<Range>>>>)
std::vector<TailIndexEstimate> hill_plot(const Range &sample, std::size_t maximumOrder = 0)
{
    using Real = std::remove_cv_t<std::ranges::range_value_t<Range>>;
    return hill_plot(std::span<const Real>(sample), maximumOrder);
}

template <std::floating_point Real>
double pickands_estimator(std::span<const Real> sample, std::size_t order)
{
    const auto sorted = detail::sortedTailSample(sample);
    if (order == 0 || order > sorted.size() / 4)
        throw std::invalid_argument(
            "Pickands estimator order should be positive with four times the order at most the sample size");
    const double upper = static_cast<double>(sorted[sorted.size() - order]);
    const double middle = static_cast<double>(sorted[sorted.size() - 2 * order]);
    const double lower = static_cast<double>(sorted[sorted.size() - 4 * order]);
    const double numerator = upper - middle;
    const double denominator = middle - lower;
    if (!(numerator > 0) || !(denominator > 0))
        throw std::invalid_argument("Pickands estimator requires distinct upper order statistics");
    return std::log(numerator / denominator) / std::log(2.0);
}

template <std::ranges::contiguous_range Range>
    requires std::floating_point<std::remove_cv_t<std::ranges::range_value_t<Range>>> &&
             (!std::same_as<std::remove_cvref_t<Range>,
                            std::span<const std::remove_cv_t<std::ranges::range_value_t<Range>>>>)
double pickands_estimator(const Range &sample, std::size_t order)
{
    using Real = std::remove_cv_t<std::ranges::range_value_t<Range>>;
    return pickands_estimator(std::span<const Real>(sample), order);
}

template <std::floating_point Real>
FitResult<GeneralizedParetoRand<Real>> fit_peaks_over_threshold(std::span<const Real> sample, Real threshold,
                                                                FitOptions options = {})
{
    if (!std::isfinite(threshold))
        throw std::invalid_argument("Peaks-over-threshold fitting requires a finite threshold");
    std::vector<Real> exceedances;
    exceedances.reserve(sample.size());
    for (Real value : sample) {
        if (!std::isfinite(value))
            throw fit_error(FitErrorCode::WrongSample, "Peaks-over-threshold fitting requires finite observations");
        if (value > threshold)
            exceedances.push_back(value - threshold);
    }
    auto result = fit<GeneralizedParetoRand<Real>>(std::span<const Real>(exceedances), std::move(options));
    result.distribution = GeneralizedParetoRand<Real>(result.distribution.GetShape(), result.distribution.GetScale(),
                                                      static_cast<double>(threshold));
    return result;
}

template <std::ranges::contiguous_range Range>
    requires std::floating_point<std::remove_cv_t<std::ranges::range_value_t<Range>>> &&
             (!std::same_as<std::remove_cvref_t<Range>,
                            std::span<const std::remove_cv_t<std::ranges::range_value_t<Range>>>>)
auto fit_peaks_over_threshold(const Range &sample, std::ranges::range_value_t<Range> threshold, FitOptions options = {})
{
    using Real = std::remove_cv_t<std::ranges::range_value_t<Range>>;
    return fit_peaks_over_threshold(std::span<const Real>(sample), static_cast<Real>(threshold), std::move(options));
}

} // namespace randlib

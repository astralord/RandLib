#pragma once

#include <randlib/fit/mle.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace randlib
{

enum class Censoring
{
    Exact,
    Left,
    Right,
    Interval
};

template <std::floating_point Real>
struct CensoredObservation
{
    Censoring censoring = Censoring::Exact;
    Real lower = 0;
    Real upper = 0;

    static CensoredObservation exact(Real value) { return {Censoring::Exact, value, value}; }
    static CensoredObservation left(Real upper) { return {Censoring::Left, -infinity(), upper}; }
    static CensoredObservation right(Real lower) { return {Censoring::Right, lower, infinity()}; }
    static CensoredObservation interval(Real lower, Real upper) { return {Censoring::Interval, lower, upper}; }

private:
    static constexpr Real infinity() { return std::numeric_limits<Real>::infinity(); }
};

namespace detail
{

template <std::floating_point Real>
void validateCensored(const CensoredObservation<Real> &observation)
{
    switch (observation.censoring) {
    case Censoring::Exact:
        if (!std::isfinite(observation.lower) || observation.lower != observation.upper)
            throw std::invalid_argument("An exact observation requires one finite value");
        break;
    case Censoring::Left:
        if (!std::isfinite(observation.upper))
            throw std::invalid_argument("A left-censored observation requires a finite upper bound");
        break;
    case Censoring::Right:
        if (!std::isfinite(observation.lower))
            throw std::invalid_argument("A right-censored observation requires a finite lower bound");
        break;
    case Censoring::Interval:
        if (!std::isfinite(observation.lower) || !std::isfinite(observation.upper) ||
            !(observation.lower < observation.upper))
            throw std::invalid_argument("An interval-censored observation requires ordered finite bounds");
        break;
    }
}

template <class Distribution>
double intervalProbability(const Distribution &distribution, typename Distribution::value_type lower,
                           typename Distribution::value_type upper)
{
    const double cdfDifference = distribution.F(upper) - distribution.F(lower);
    const double survivalDifference = distribution.S(lower) - distribution.S(upper);
    return std::max(cdfDifference, survivalDifference);
}

} // namespace detail

template <class Distribution>
    requires std::floating_point<typename Distribution::value_type> && detail::FittableModel<Distribution>
FitResult<Distribution>
fit_censored(std::span<const CensoredObservation<typename Distribution::value_type>> observations,
             FitOptions options = {})
{
    using Real = typename Distribution::value_type;
    if (observations.empty())
        throw fit_error(FitErrorCode::TooFewElements, "Censored fitting requires a non-empty sample");
    if (options.method != Method::MLE)
        throw fit_error(FitErrorCode::NotApplicable, "Censored fitting supports maximum likelihood only");

    std::vector<Real> initializationSample;
    initializationSample.reserve(observations.size());
    for (const auto &observation : observations) {
        detail::validateCensored(observation);
        switch (observation.censoring) {
        case Censoring::Exact:
        case Censoring::Right:
            initializationSample.push_back(observation.lower);
            break;
        case Censoring::Left:
            initializationSample.push_back(observation.upper);
            break;
        case Censoring::Interval:
            initializationSample.push_back(observation.lower + (observation.upper - observation.lower) / 2);
            break;
        }
    }

    const auto logLikelihood = [observations](const Distribution &distribution) {
        long double result = 0;
        for (const auto &observation : observations) {
            switch (observation.censoring) {
            case Censoring::Exact:
                result += distribution.logf(observation.lower);
                break;
            case Censoring::Left: {
                const double probability = distribution.F(observation.upper);
                if (!(probability > 0))
                    return -std::numeric_limits<double>::infinity();
                result += std::log(probability);
                break;
            }
            case Censoring::Right: {
                const double probability = distribution.S(observation.lower);
                if (!(probability > 0))
                    return -std::numeric_limits<double>::infinity();
                result += std::log(probability);
                break;
            }
            case Censoring::Interval: {
                const double probability =
                    detail::intervalProbability(distribution, observation.lower, observation.upper);
                if (!(probability > 0))
                    return -std::numeric_limits<double>::infinity();
                result += std::log(probability);
                break;
            }
            }
        }
        return static_cast<double>(result);
    };
    return detail::fitLikelihood<Distribution>(std::span<const Real>(initializationSample), std::move(options),
                                               logLikelihood, true);
}

template <class Distribution, std::ranges::contiguous_range Range>
    requires std::floating_point<typename Distribution::value_type> && detail::FittableModel<Distribution> &&
             std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>,
                          CensoredObservation<typename Distribution::value_type>> &&
             (!std::same_as<std::remove_cvref_t<Range>,
                            std::span<const CensoredObservation<typename Distribution::value_type>>>)
FitResult<Distribution> fit_censored(const Range &observations, FitOptions options = {})
{
    using Observation = CensoredObservation<typename Distribution::value_type>;
    return fit_censored<Distribution>(std::span<const Observation>(observations), std::move(options));
}

template <class Distribution>
    requires std::floating_point<typename Distribution::value_type> && detail::FittableModel<Distribution>
FitResult<Distribution> fit_truncated(std::span<const typename Distribution::value_type> sample,
                                      typename Distribution::value_type lower, typename Distribution::value_type upper,
                                      FitOptions options = {})
{
    using Real = typename Distribution::value_type;
    if (options.method != Method::MLE)
        throw fit_error(FitErrorCode::NotApplicable, "Truncated fitting supports maximum likelihood only");
    if (std::isnan(lower) || std::isnan(upper) || !(lower < upper))
        throw std::invalid_argument("Truncation bounds should be ordered and not NaN");
    if (std::ranges::any_of(
            sample, [lower, upper](Real value) { return !std::isfinite(value) || value < lower || value > upper; }))
        throw fit_error(FitErrorCode::WrongSample, "Truncated observations should lie within the truncation bounds");

    const auto logLikelihood = [sample, lower, upper](const Distribution &distribution) {
        const double lowerCdf = std::isinf(lower) ? 0 : distribution.F(lower);
        const double upperCdf = std::isinf(upper) ? 1 : distribution.F(upper);
        const double lowerSurvival = std::isinf(lower) ? 1 : distribution.S(lower);
        const double upperSurvival = std::isinf(upper) ? 0 : distribution.S(upper);
        const double probability = std::max(upperCdf - lowerCdf, lowerSurvival - upperSurvival);
        if (!(probability > 0))
            return -std::numeric_limits<double>::infinity();
        return distribution.LogLikelihoodFunction(sample) - sample.size() * std::log(probability);
    };
    return detail::fitLikelihood<Distribution>(sample, std::move(options), logLikelihood, true);
}

template <class Distribution, std::ranges::contiguous_range Range>
    requires std::floating_point<typename Distribution::value_type> && detail::FittableModel<Distribution> &&
             std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, typename Distribution::value_type> &&
             (!std::same_as<std::remove_cvref_t<Range>, std::span<const typename Distribution::value_type>>)
FitResult<Distribution> fit_truncated(const Range &sample, typename Distribution::value_type lower,
                                      typename Distribution::value_type upper, FitOptions options = {})
{
    using Real = typename Distribution::value_type;
    return fit_truncated<Distribution>(std::span<const Real>(sample), lower, upper, std::move(options));
}

} // namespace randlib

#pragma once

#include <randlib/fit/mle.hpp>

#include <cmath>
#include <concepts>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>

namespace randlib
{

template <class... Distributions>
struct BestFitResult
{
    std::variant<FitResult<Distributions>...> fit;
    std::size_t index = 0;
    Criterion criterion = Criterion::AIC;
    double value = std::numeric_limits<double>::infinity();
};

template <class First, class... Rest>
    requires(detail::FittableModel<First> && ... && detail::FittableModel<Rest>) &&
            (std::same_as<typename First::value_type, typename Rest::value_type> && ...)
BestFitResult<First, Rest...> fit_best(std::span<const typename First::value_type> sample, Criterion criterion)
{
    using Result = BestFitResult<First, Rest...>;
    using Variant = decltype(Result::fit);
    std::optional<Variant> best;
    double bestValue = std::numeric_limits<double>::infinity();
    std::size_t bestIndex = 0;
    std::size_t index = 0;
    auto consider = [&]<class Distribution>() {
        try {
            auto candidate = fit<Distribution>(sample);
            const double value = criterion == Criterion::AIC ? candidate.aic : candidate.bic;
            if (candidate.converged && std::isfinite(value) && value < bestValue) {
                best.emplace(std::move(candidate));
                bestValue = value;
                bestIndex = index;
            }
        } catch (const std::exception &) {
        }
        ++index;
    };
    consider.template operator()<First>();
    (consider.template operator()<Rest>(), ...);
    if (!best)
        throw fit_error(FitErrorCode::Undefined, "No candidate distribution produced a finite converged fit");
    return {std::move(*best), bestIndex, criterion, bestValue};
}

template <class First, class... Rest, std::ranges::contiguous_range Range>
    requires std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, typename First::value_type> &&
             (detail::FittableModel<First> && ... && detail::FittableModel<Rest>) &&
             (std::same_as<typename First::value_type, typename Rest::value_type> && ...)
BestFitResult<First, Rest...> fit_best(const Range &sample, Criterion criterion)
{
    return fit_best<First, Rest...>(std::span<const typename First::value_type>(sample), criterion);
}

} // namespace randlib

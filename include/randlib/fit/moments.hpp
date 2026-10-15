#pragma once

#include <randlib/fit/mle.hpp>

#include <ranges>
#include <span>
#include <type_traits>

namespace randlib
{

template <class Distribution>
    requires detail::FittableModel<Distribution>
FitResult<Distribution> fit_moments(std::span<const typename Distribution::value_type> sample)
{
    return fit<Distribution>(sample, Method::MethodOfMoments);
}

template <class Distribution, std::ranges::contiguous_range Range>
    requires std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, typename Distribution::value_type> &&
             detail::FittableModel<Distribution>
FitResult<Distribution> fit_moments(const Range &sample)
{
    return fit_moments<Distribution>(std::span<const typename Distribution::value_type>(sample));
}

} // namespace randlib

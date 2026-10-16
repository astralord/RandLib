#pragma once

#include <randlib/core/errors.hpp>
#include <randlib/bivariate/normal_inverse_gamma.hpp>
#include <randlib/continuous/beta.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/inverse_gamma.hpp>
#include <randlib/continuous/log_normal.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/pareto.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/weibull.hpp>
#include <randlib/discrete/binomial.hpp>
#include <randlib/discrete/negative_binomial.hpp>
#include <randlib/discrete/poisson.hpp>

#include <concepts>
#include <ranges>
#include <span>
#include <type_traits>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif

namespace randlib
{

namespace detail
{

template <class Distribution, class Prior>
struct BayesModel;

template <class Distribution, class Prior>
struct ConditionalBayesModel;

template <class Distribution, class Prior>
concept ConjugateModel = requires(std::span<const typename Distribution::value_type> sample, const Prior &prior) {
    BayesModel<Distribution, Prior>::posterior(sample, prior);
};

template <class Distribution, class Prior, class Fixed>
concept ConditionalConjugateModel =
    requires(std::span<const typename Distribution::value_type> sample, const Prior &prior, Fixed fixed) {
        ConditionalBayesModel<Distribution, Prior>::posterior(sample, prior, fixed);
    };

template <std::floating_point Real>
struct BayesModel<NormalRand<Real>, NormalInverseGammaRand<Real>>
{
    static NormalInverseGammaRand<Real> posterior(std::span<const Real> sample,
                                                  const NormalInverseGammaRand<Real> &prior)
    {
        NormalRand<Real> distribution;
        return distribution.FitBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct BayesModel<LogNormalRand<Real>, NormalInverseGammaRand<Real>>
{
    static NormalInverseGammaRand<Real> posterior(std::span<const Real> sample,
                                                  const NormalInverseGammaRand<Real> &prior)
    {
        LogNormalRand<Real> distribution;
        return distribution.FitBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct BayesModel<ExponentialRand<Real>, GammaRand<Real>>
{
    static GammaRand<Real> posterior(std::span<const Real> sample, const GammaRand<Real> &prior)
    {
        ExponentialRand<Real> distribution;
        return distribution.FitRateBayes(sample, prior);
    }
};

template <std::integral Int>
struct BayesModel<PoissonRand<Int>, GammaRand<double>>
{
    static GammaRand<double> posterior(std::span<const Int> sample, const GammaRand<double> &prior)
    {
        PoissonRand<Int> distribution;
        return distribution.FitBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct ConditionalBayesModel<NormalRand<Real>, NormalRand<Real>>
{
    static NormalRand<Real> posterior(std::span<const Real> sample, const NormalRand<Real> &prior, double variance)
    {
        NormalRand<Real> distribution(0, variance);
        return distribution.FitLocationBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct ConditionalBayesModel<NormalRand<Real>, InverseGammaRand<Real>>
{
    static InverseGammaRand<Real> posterior(std::span<const Real> sample, const InverseGammaRand<Real> &prior,
                                            double location)
    {
        NormalRand<Real> distribution(location, 1);
        return distribution.FitVarianceBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct ConditionalBayesModel<LogNormalRand<Real>, NormalRand<Real>>
{
    static NormalRand<Real> posterior(std::span<const Real> sample, const NormalRand<Real> &prior, double variance)
    {
        LogNormalRand<Real> distribution(0, variance);
        return distribution.FitLocationBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct ConditionalBayesModel<LogNormalRand<Real>, InverseGammaRand<Real>>
{
    static InverseGammaRand<Real> posterior(std::span<const Real> sample, const InverseGammaRand<Real> &prior,
                                            double location)
    {
        LogNormalRand<Real> distribution(location, 1);
        return distribution.FitScaleBayes(sample, prior);
    }
};

template <std::integral Int>
struct ConditionalBayesModel<BinomialRand<Int>, BetaRand<double>>
{
    static BetaRand<double> posterior(std::span<const Int> sample, const BetaRand<double> &prior, int number)
    {
        BinomialRand<Int> distribution(number);
        return distribution.FitProbabilityBayes(sample, prior);
    }
};

template <std::integral Int, class Shape>
struct ConditionalBayesModel<NegativeBinomialRand<Int, Shape>, BetaRand<double>>
{
    static BetaRand<double> posterior(std::span<const Int> sample, const BetaRand<double> &prior, Shape number)
    {
        NegativeBinomialRand<Int, Shape> distribution(number);
        return distribution.FitProbabilityBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct ConditionalBayesModel<WeibullRand<Real>, InverseGammaRand<Real>>
{
    static InverseGammaRand<Real> posterior(std::span<const Real> sample, const InverseGammaRand<Real> &prior,
                                            double shape)
    {
        WeibullRand<Real> distribution(1, shape);
        return distribution.FitScaleBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct ConditionalBayesModel<ParetoRand<Real>, GammaRand<Real>>
{
    static GammaRand<Real> posterior(std::span<const Real> sample, const GammaRand<Real> &prior, double scale)
    {
        ParetoRand<Real> distribution(1, scale);
        return distribution.FitShapeBayes(sample, prior);
    }
};

template <std::floating_point Real>
struct ConditionalBayesModel<UniformRand<Real>, ParetoRand<Real>>
{
    static ParetoRand<Real> posterior(std::span<const Real> sample, const ParetoRand<Real> &prior, double minimum)
    {
        UniformRand<Real> distribution(minimum, minimum + 1);
        return distribution.FitMaximumBayes(sample, prior);
    }
};

} // namespace detail

template <class Distribution, class Prior>
    requires detail::ConjugateModel<Distribution, Prior>
auto fit_bayes(std::span<const typename Distribution::value_type> sample, const Prior &prior)
{
    if (sample.empty())
        throw fit_error(FitErrorCode::TooFewElements, "Bayesian fitting requires a non-empty sample");
    return detail::BayesModel<Distribution, Prior>::posterior(sample, prior);
}

template <class Distribution, std::ranges::contiguous_range Range, class Prior>
    requires std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, typename Distribution::value_type> &&
             detail::ConjugateModel<Distribution, Prior>
auto fit_bayes(const Range &sample, const Prior &prior)
{
    return fit_bayes<Distribution>(std::span<const typename Distribution::value_type>(sample), prior);
}

template <class Distribution, std::ranges::contiguous_range Range, class Prior, class Fixed>
    requires std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, typename Distribution::value_type> &&
             detail::ConditionalConjugateModel<Distribution, Prior, Fixed>
auto fit_bayes(const Range &sample, const Prior &prior, Fixed fixed)
{
    if (std::ranges::empty(sample))
        throw fit_error(FitErrorCode::TooFewElements, "Bayesian fitting requires a non-empty sample");
    return detail::ConditionalBayesModel<Distribution, Prior>::posterior(
        std::span<const typename Distribution::value_type>(sample), prior, fixed);
}

template <std::ranges::contiguous_range Range, std::floating_point Real>
    requires std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, Real>
auto fit_bayes(const Range &sample, const NormalInverseGammaRand<Real> &prior)
{
    return fit_bayes<NormalRand<Real>>(sample, prior);
}

template <std::ranges::contiguous_range Range, std::floating_point Real>
    requires std::same_as<std::remove_cv_t<std::ranges::range_value_t<Range>>, Real>
auto fit_bayes(const Range &sample, const GammaRand<Real> &prior)
{
    return fit_bayes<ExponentialRand<Real>>(sample, prior);
}

template <std::ranges::contiguous_range Range>
    requires std::integral<std::remove_cv_t<std::ranges::range_value_t<Range>>>
auto fit_bayes(const Range &sample, const GammaRand<double> &prior)
{
    using Int = std::remove_cv_t<std::ranges::range_value_t<Range>>;
    return fit_bayes<PoissonRand<Int>>(sample, prior);
}

} // namespace randlib

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif

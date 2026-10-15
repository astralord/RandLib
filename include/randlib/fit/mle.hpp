#pragma once

#include <randlib/continuous/beta.hpp>
#include <randlib/continuous/beta_prime.hpp>
#include <randlib/continuous/degenerate.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/laplace.hpp>
#include <randlib/continuous/levy.hpp>
#include <randlib/continuous/log_normal.hpp>
#include <randlib/continuous/logistic.hpp>
#include <randlib/continuous/nakagami.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/pareto.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/weibull.hpp>
#include <randlib/discrete/binomial.hpp>
#include <randlib/discrete/negative_binomial.hpp>
#include <randlib/discrete/poisson.hpp>
#include <randlib/discrete/uniform_discrete.hpp>
#include <randlib/fit/result.hpp>
#include <randlib/math/optimize.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <initializer_list>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace randlib
{

namespace detail
{

struct FitParameter
{
    std::string_view name;
    double value;
    std::optional<double> lower;
    std::optional<double> upper;
    bool lowerInclusive = false;
    bool upperInclusive = false;
    bool integral = false;
};

using FixedParameterMap = std::map<std::string, double, std::less<>>;

template <class Distribution>
struct FitModel;

template <class Distribution>
concept FittableModel = requires(Distribution distribution, std::span<const typename Distribution::value_type> sample,
                                 std::span<const double> parameters, Method method, const FixedParameterMap &fixed) {
    { FitModel<Distribution>::parameters(distribution) } -> std::same_as<std::vector<FitParameter>>;
    { FitModel<Distribution>::make(parameters) } -> std::same_as<Distribution>;
    { FitModel<Distribution>::initial(sample, method, fixed) } -> std::same_as<Distribution>;
    { FitModel<Distribution>::hasClosedForm(method) } -> std::same_as<bool>;
};

inline double toUnconstrained(double value, const FitParameter &parameter)
{
    if (parameter.lower && parameter.upper) {
        const double fraction = (value - *parameter.lower) / (*parameter.upper - value);
        return std::log(fraction);
    }
    if (parameter.lower)
        return std::log(value - *parameter.lower);
    if (parameter.upper)
        return std::log(*parameter.upper - value);
    return value;
}

inline double fromUnconstrained(double value, const FitParameter &parameter)
{
    if (parameter.lower && parameter.upper) {
        const double width = *parameter.upper - *parameter.lower;
        if (value >= 0) {
            const double inverse = std::exp(-value);
            return *parameter.lower + width / (1 + inverse);
        }
        const double exponential = std::exp(value);
        return *parameter.lower + width * exponential / (1 + exponential);
    }
    if (parameter.lower)
        return *parameter.lower + std::exp(value);
    if (parameter.upper)
        return *parameter.upper - std::exp(value);
    return value;
}

inline bool withinBounds(double value, const FitParameter &parameter)
{
    if (!std::isfinite(value))
        return false;
    if (parameter.lower && (value < *parameter.lower || (!parameter.lowerInclusive && value == *parameter.lower)))
        return false;
    if (parameter.upper && (value > *parameter.upper || (!parameter.upperInclusive && value == *parameter.upper)))
        return false;
    if (parameter.integral && value != std::round(value))
        return false;
    return true;
}

inline void requireMethod(Method method, std::string_view distribution,
                          std::initializer_list<Method> supported = {Method::MLE, Method::MethodOfMoments})
{
    if (std::ranges::find(supported, method) == supported.end())
        throw fit_error(FitErrorCode::NotApplicable,
                        std::string(distribution) + " distribution: requested fitting method is not implemented");
}

inline std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> matrix)
{
    const std::size_t n = matrix.size();
    std::vector<std::vector<double>> result(n, std::vector<double>(n));
    for (std::size_t i = 0; i < n; ++i)
        result[i][i] = 1;

    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < n; ++row)
            if (std::fabs(matrix[row][column]) > std::fabs(matrix[pivot][column]))
                pivot = row;
        if (!std::isfinite(matrix[pivot][column]) ||
            std::fabs(matrix[pivot][column]) <= std::numeric_limits<double>::epsilon())
            return {};
        std::swap(matrix[column], matrix[pivot]);
        std::swap(result[column], result[pivot]);

        const double diagonal = matrix[column][column];
        for (std::size_t j = 0; j < n; ++j) {
            matrix[column][j] /= diagonal;
            result[column][j] /= diagonal;
        }
        for (std::size_t row = 0; row < n; ++row) {
            if (row == column)
                continue;
            const double factor = matrix[row][column];
            for (std::size_t j = 0; j < n; ++j) {
                matrix[row][j] -= factor * matrix[column][j];
                result[row][j] -= factor * result[column][j];
            }
        }
    }
    return result;
}

template <class Distribution>
double negativeLogLikelihood(const std::vector<double> &parameters,
                             std::span<const typename Distribution::value_type> sample)
{
    try {
        Distribution distribution = FitModel<Distribution>::make(parameters);
        const double value = -distribution.LogLikelihoodFunction(sample);
        return std::isfinite(value) ? value : std::numeric_limits<double>::infinity();
    } catch (const std::exception &) {
        return std::numeric_limits<double>::infinity();
    }
}

template <class Distribution>
void addUncertainty(FitResult<Distribution> &result, std::span<const typename Distribution::value_type> sample)
{
    std::vector<std::size_t> freeIndices;
    std::vector<double> freeParameters;
    for (std::size_t i = 0; i < result.parameters.size(); ++i) {
        if (!result.fixed[i]) {
            freeIndices.push_back(i);
            freeParameters.push_back(result.parameters[i]);
        }
    }

    const std::size_t dimension = result.parameters.size();
    result.covariance.assign(dimension, std::vector<double>(dimension));
    result.std_errors.assign(dimension, 0);
    if (freeIndices.empty())
        return;

    auto objective = [&](std::span<const double> free) {
        std::vector<double> parameters = result.parameters;
        for (std::size_t i = 0; i < freeIndices.size(); ++i)
            parameters[freeIndices[i]] = free[i];
        return negativeLogLikelihood<Distribution>(parameters, sample);
    };
    auto information = math::numericalHessian(objective, freeParameters);
    auto covariance = inverse(std::move(information));
    if (covariance.empty()) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t index : freeIndices) {
            result.std_errors[index] = nan;
            for (std::size_t other : freeIndices)
                result.covariance[index][other] = nan;
        }
        return;
    }
    for (std::size_t i = 0; i < freeIndices.size(); ++i) {
        const std::size_t row = freeIndices[i];
        result.std_errors[row] =
            covariance[i][i] >= 0 ? std::sqrt(covariance[i][i]) : std::numeric_limits<double>::quiet_NaN();
        for (std::size_t j = 0; j < freeIndices.size(); ++j)
            result.covariance[row][freeIndices[j]] = covariance[i][j];
    }
}

template <std::floating_point Real>
struct FitModel<NormalRand<Real>>
{
    using Distribution = NormalRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"location", distribution.GetLocation(), std::nullopt, std::nullopt},
                {"variance", static_cast<double>(distribution.Variance()), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        Distribution distribution;
        if (method == Method::MLE || method == Method::MethodOfMoments)
            distribution.Fit(sample);
        else if (method == Method::UMVU)
            distribution.Fit(sample, true);
        else
            throw fit_error(FitErrorCode::NotApplicable, "Normal distribution: quantile fitting is not implemented");
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<LogNormalRand<Real>>
{
    using Distribution = LogNormalRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        const double scale = distribution.GetScale();
        return {{"location", distribution.GetLocation(), std::nullopt, std::nullopt},
                {"variance", scale * scale, 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Log-normal", {Method::MLE});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<ExponentialRand<Real>>
{
    using Distribution = ExponentialRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"rate", distribution.GetRate(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Exponential", {Method::MLE, Method::UMVU});
        Distribution distribution;
        distribution.FitRate(sample, method == Method::UMVU);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<GammaRand<Real>>
{
    using Distribution = GammaRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), 0.0, std::nullopt},
                {"rate", distribution.GetRate(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Gamma", {Method::MLE});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<BetaRand<Real>>
{
    using Distribution = BetaRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"alpha", distribution.GetAlpha(), 0.0, std::nullopt},
                {"beta", distribution.GetBeta(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Beta", {Method::MLE});
        Distribution distribution;
        distribution.FitShapes(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<ArcsineRand<Real>>
{
    using Distribution = ArcsineRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Arcsine", {Method::MLE});
        Distribution distribution;
        distribution.FitShape(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<BetaPrimeRand<Real>>
{
    using Distribution = BetaPrimeRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"alpha", distribution.GetAlpha(), 0.0, std::nullopt},
                {"beta", distribution.GetBeta(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Beta prime", {Method::MLE});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<ParetoRand<Real>>
{
    using Distribution = ParetoRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), 0.0, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Pareto", {Method::MLE, Method::UMVU});
        Distribution distribution;
        distribution.Fit(sample, method == Method::UMVU);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<WeibullRand<Real>>
{
    using Distribution = WeibullRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"scale", distribution.GetScale(), 0.0, std::nullopt},
                {"shape", distribution.GetShape(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Weibull", {Method::MLE});
        if (std::ranges::any_of(sample, [](Real value) { return !(value > 0); }))
            throw fit_error(FitErrorCode::WrongSample, "Weibull distribution: sample values should be positive");
        const long double sum = std::accumulate(sample.begin(), sample.end(), 0.0L);
        return Distribution(static_cast<double>(sum / sample.size()), 1.0);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<LaplaceRand<Real>>
{
    using Distribution = LaplaceRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"location", distribution.GetShift(), std::nullopt, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Laplace", {Method::MLE});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<LogisticRand<Real>>
{
    using Distribution = LogisticRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"location", distribution.GetLocation(), std::nullopt, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Logistic", {Method::MLE});
        const long double sum = std::accumulate(sample.begin(), sample.end(), 0.0L);
        const double location = static_cast<double>(sum / sample.size());
        long double squaredDeviation = 0;
        for (Real value : sample) {
            const long double difference = value - location;
            squaredDeviation += difference * difference;
        }
        const double scale = std::sqrt(static_cast<double>(3 * squaredDeviation / sample.size())) / std::numbers::pi;
        if (!(scale > 0))
            throw fit_error(FitErrorCode::WrongSample, "Logistic distribution: sample variance should be positive");
        return Distribution(location, scale);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<LevyRand<Real>>
{
    using Distribution = LevyRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"location", distribution.GetLocation(), std::nullopt, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Levy", {Method::MLE});
        const auto location = fixed.find("location");
        if (location == fixed.end())
            throw fit_error(FitErrorCode::NotApplicable,
                            "Levy distribution: fitting currently requires a fixed location");
        Distribution distribution(location->second, 1.0);
        distribution.FitScale(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<RayleighRand<Real>>
{
    using Distribution = RayleighRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Rayleigh", {Method::MLE, Method::UMVU});
        Distribution distribution;
        distribution.Fit(sample, method == Method::UMVU);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<UniformRand<Real>>
{
    using Distribution = UniformRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"minimum", distribution.MinValue(), std::nullopt, distribution.MaxValue()},
                {"maximum", distribution.MaxValue(), distribution.MinValue(), std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Uniform", {Method::MLE, Method::UMVU});
        Distribution distribution;
        distribution.Fit(sample, method == Method::UMVU);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <>
struct FitModel<DegenerateRand>
{
    using Distribution = DegenerateRand;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"value", distribution.GetValue(), std::nullopt, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0]); }

    static Distribution initial(std::span<const double> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Degenerate", {Method::MLE});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::integral Int>
struct FitModel<PoissonRand<Int>>
{
    using Distribution = PoissonRand<Int>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"rate", distribution.GetRate(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0]); }

    static Distribution initial(std::span<const Int> sample, Method method, const FixedParameterMap &)
    {
        if (method == Method::Quantile)
            throw fit_error(FitErrorCode::NotApplicable, "Poisson distribution: quantile fitting is not implemented");
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::integral Int>
struct FitModel<BinomialRand<Int>>
{
    using Distribution = BinomialRand<Int>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"number", static_cast<double>(distribution.GetNumber()), 1.0, std::nullopt, true, false, true},
                {"probability", distribution.GetProbability(), 0.0, 1.0, true, true}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(static_cast<int>(parameters[0]), parameters[1]);
    }

    static Distribution initial(std::span<const Int> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Binomial", {Method::MLE});
        const auto number = fixed.find("number");
        if (number == fixed.end())
            throw fit_error(FitErrorCode::NotApplicable,
                            "Binomial distribution: fitting currently requires a fixed number");
        if (number->second < 1 || number->second != std::round(number->second))
            throw std::invalid_argument("Fixed binomial number should be a positive integer");
        Distribution distribution(static_cast<int>(number->second));
        distribution.FitProbability(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::integral Int, std::floating_point Shape>
struct FitModel<NegativeBinomialRand<Int, Shape>>
{
    using Distribution = NegativeBinomialRand<Int, Shape>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"number", distribution.GetNumber(), 0.0, std::nullopt},
                {"probability", distribution.GetProbability(), 0.0, 1.0}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(static_cast<Shape>(parameters[0]), parameters[1]);
    }

    static Distribution initial(std::span<const Int> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Negative binomial", {Method::MLE});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::integral Int>
struct FitModel<UniformDiscreteRand<Int>>
{
    using Distribution = UniformDiscreteRand<Int>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"minimum", static_cast<double>(distribution.MinValue()), std::nullopt,
                 static_cast<double>(distribution.MaxValue()), false, false, true},
                {"maximum", static_cast<double>(distribution.MaxValue()), static_cast<double>(distribution.MinValue()),
                 std::nullopt, false, false, true}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(static_cast<Int>(parameters[0]), static_cast<Int>(parameters[1]));
    }

    static Distribution initial(std::span<const Int> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Uniform discrete", {Method::MLE});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

} // namespace detail

template <class Distribution>
    requires detail::FittableModel<Distribution>
FitResult<Distribution> fit(std::span<const typename Distribution::value_type> sample, FitOptions options = {})
{
    if (sample.empty())
        throw fit_error(FitErrorCode::TooFewElements, "Fitting requires a non-empty sample");
    if (!(options.tolerance > 0) || options.max_iterations <= 0)
        throw std::invalid_argument("Fitting tolerance and iteration limit should be positive");
    if (options.method != Method::MLE && !options.fixed.empty())
        throw fit_error(FitErrorCode::NotApplicable, "Fixed parameters are currently supported only for MLE");

    Distribution initial = detail::FitModel<Distribution>::initial(sample, options.method, options.fixed);
    auto specifications = detail::FitModel<Distribution>::parameters(initial);
    std::vector<double> parameters;
    parameters.reserve(specifications.size());
    std::vector<bool> fixed(specifications.size());
    for (std::size_t i = 0; i < specifications.size(); ++i) {
        auto found = options.fixed.find(specifications[i].name);
        if (found != options.fixed.end()) {
            if (!detail::withinBounds(found->second, specifications[i]))
                throw std::invalid_argument("Fixed parameter is outside its bounds: " +
                                            std::string(specifications[i].name));
            specifications[i].value = found->second;
            fixed[i] = true;
        }
        parameters.push_back(specifications[i].value);
    }
    for (const auto &[name, value] : options.fixed) {
        (void)value;
        if (std::ranges::none_of(specifications,
                                 [&name](const detail::FitParameter &parameter) { return parameter.name == name; }))
            throw std::invalid_argument("Unknown fixed parameter: " + name);
    }

    std::vector<std::size_t> freeIndices;
    std::vector<double> unconstrained;
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        if (!fixed[i]) {
            if (specifications[i].integral &&
                (!detail::FitModel<Distribution>::hasClosedForm(options.method) || !options.fixed.empty()))
                throw fit_error(FitErrorCode::NotApplicable,
                                "Free integer parameters are not supported by continuous optimization");
            freeIndices.push_back(i);
            unconstrained.push_back(detail::toUnconstrained(parameters[i], specifications[i]));
        }
    }

    math::MinimizeResult<double> optimization;
    optimization.x = unconstrained;
    optimization.converged = true;
    if (!freeIndices.empty() &&
        (!detail::FitModel<Distribution>::hasClosedForm(options.method) || !options.fixed.empty())) {
        auto objective = [&](std::span<const double> free) {
            std::vector<double> candidate = parameters;
            for (std::size_t i = 0; i < free.size(); ++i)
                candidate[freeIndices[i]] = detail::fromUnconstrained(free[i], specifications[freeIndices[i]]);
            return detail::negativeLogLikelihood<Distribution>(candidate, sample);
        };
        optimization = math::bfgs(objective, unconstrained, options.tolerance, options.max_iterations);
        if (!optimization.converged)
            optimization =
                math::nelderMead(objective, std::move(optimization.x), options.tolerance, options.max_iterations);
        for (std::size_t i = 0; i < optimization.x.size(); ++i)
            parameters[freeIndices[i]] = detail::fromUnconstrained(optimization.x[i], specifications[freeIndices[i]]);
    }

    Distribution distribution = detail::FitModel<Distribution>::make(parameters);
    const double logLikelihood = distribution.LogLikelihoodFunction(sample);
    const std::size_t fittedParameters = freeIndices.size();
    FitResult<Distribution> result(std::move(distribution));
    result.log_likelihood = logLikelihood;
    result.aic = 2 * fittedParameters - 2 * logLikelihood;
    result.bic = std::log(static_cast<double>(sample.size())) * fittedParameters - 2 * logLikelihood;
    result.parameters = std::move(parameters);
    result.fixed = std::move(fixed);
    result.converged = optimization.converged;
    result.iterations = optimization.iterations;
    result.method = options.method;
    for (const auto &parameter : specifications)
        result.parameter_names.emplace_back(parameter.name);
    detail::addUncertainty(result, sample);
    return result;
}

template <class Distribution>
    requires detail::FittableModel<Distribution>
FitResult<Distribution> fit(std::span<const typename Distribution::value_type> sample, Method method)
{
    FitOptions options;
    options.method = method;
    return fit<Distribution>(sample, std::move(options));
}

} // namespace randlib

#pragma once

#include <randlib/continuous/beta.hpp>
#include <randlib/continuous/beta_prime.hpp>
#include <randlib/continuous/cauchy.hpp>
#include <randlib/continuous/degenerate.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/continuous/fisher_f.hpp>
#include <randlib/continuous/frechet.hpp>
#include <randlib/continuous/generalized_extreme_value.hpp>
#include <randlib/continuous/generalized_pareto.hpp>
#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/gumbel.hpp>
#include <randlib/continuous/inverse_gamma.hpp>
#include <randlib/continuous/inverse_gaussian.hpp>
#include <randlib/continuous/laplace.hpp>
#include <randlib/continuous/levy.hpp>
#include <randlib/continuous/log_normal.hpp>
#include <randlib/continuous/logistic.hpp>
#include <randlib/continuous/log_stable.hpp>
#include <randlib/continuous/nakagami.hpp>
#include <randlib/continuous/noncentral_chi_squared.hpp>
#include <randlib/continuous/noncentral_t.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/pareto.hpp>
#include <randlib/continuous/student_t.hpp>
#include <randlib/continuous/stable.hpp>
#include <randlib/continuous/stable_count.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/weibull.hpp>
#include <randlib/discrete/bernoulli.hpp>
#include <randlib/discrete/binomial.hpp>
#include <randlib/discrete/geometric.hpp>
#include <randlib/discrete/negative_binomial.hpp>
#include <randlib/discrete/poisson.hpp>
#include <randlib/discrete/uniform_discrete.hpp>
#include <randlib/fit/result.hpp>
#include <randlib/math/optimize.hpp>

#include <algorithm>
#include <array>
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

template <class Distribution>
bool requiresOptimization(Method method, const FixedParameterMap &fixed)
{
    if constexpr (requires { FitModel<Distribution>::requiresOptimization(method, fixed); })
        return FitModel<Distribution>::requiresOptimization(method, fixed);
    return !FitModel<Distribution>::hasClosedForm(method) || !fixed.empty();
}

template <class T>
DoublePair sampleMeanAndVariance(std::span<const T> sample)
{
    const long double sum = std::accumulate(sample.begin(), sample.end(), 0.0L);
    const long double mean = sum / sample.size();
    long double squaredDeviation = 0;
    for (T value : sample) {
        const long double difference = value - mean;
        squaredDeviation += difference * difference;
    }
    return {static_cast<double>(mean), static_cast<double>(squaredDeviation / sample.size())};
}

template <class T>
double sampleQuantile(std::span<const T> sample, double probability)
{
    std::vector<T> sorted(sample.begin(), sample.end());
    std::ranges::sort(sorted);
    const double index = probability * (sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(index);
    const std::size_t upper = std::min(lower + 1, sorted.size() - 1);
    const double fraction = index - lower;
    return static_cast<double>(sorted[lower]) * (1 - fraction) + static_cast<double>(sorted[upper]) * fraction;
}

template <std::size_t Rows, std::size_t Columns>
double bilinearInterpolation(const std::array<double, Rows> &rowGrid, const std::array<double, Columns> &columnGrid,
                             const std::array<std::array<double, Columns>, Rows> &values, double row, double column)
{
    const auto bracket = [](const auto &grid, double value) {
        const auto upper = std::ranges::upper_bound(grid, value);
        if (upper == grid.begin())
            return std::pair<std::size_t, double>{0, 0};
        if (upper == grid.end())
            return std::pair<std::size_t, double>{grid.size() - 2, 1};
        const std::size_t lowerIndex = static_cast<std::size_t>(upper - grid.begin() - 1);
        const double fraction = (value - grid[lowerIndex]) / (grid[lowerIndex + 1] - grid[lowerIndex]);
        return std::pair<std::size_t, double>{lowerIndex, fraction};
    };
    const auto [i, rowFraction] = bracket(rowGrid, row);
    const auto [j, columnFraction] = bracket(columnGrid, column);
    const double lower = values[i][j] * (1 - columnFraction) + values[i][j + 1] * columnFraction;
    const double upper = values[i + 1][j] * (1 - columnFraction) + values[i + 1][j + 1] * columnFraction;
    return lower * (1 - rowFraction) + upper * rowFraction;
}

inline std::optional<int> fixedPositiveInteger(const FixedParameterMap &fixed, std::string_view name,
                                               std::string_view distribution)
{
    const auto found = fixed.find(name);
    if (found == fixed.end())
        return std::nullopt;
    if (found->second < 1 || found->second > std::numeric_limits<int>::max() ||
        found->second != std::round(found->second))
        throw std::invalid_argument("Fixed " + std::string(distribution) + " " + std::string(name) +
                                    " should be a positive integer");
    return static_cast<int>(found->second);
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

template <class Distribution, class LogLikelihood>
double negativeLogLikelihood(std::span<const double> parameters, const LogLikelihood &logLikelihood)
{
    try {
        Distribution distribution = FitModel<Distribution>::make(parameters);
        const double value = -logLikelihood(distribution);
        return std::isfinite(value) ? value : std::numeric_limits<double>::infinity();
    } catch (const std::exception &) {
        return std::numeric_limits<double>::infinity();
    }
}

template <class Distribution, class Objective>
void addUncertainty(FitResult<Distribution> &result, std::span<const FitParameter> specifications,
                    const Objective &naturalObjective)
{
    const std::size_t dimension = result.parameters.size();
    result.covariance.assign(dimension, std::vector<double>(dimension));
    result.std_errors.assign(dimension, 0);
    std::vector<std::size_t> freeIndices;
    std::vector<double> freeParameters;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < result.parameters.size(); ++i) {
        if (!result.fixed[i]) {
            if (specifications[i].integral) {
                result.std_errors[i] = nan;
                std::ranges::fill(result.covariance[i], nan);
                for (auto &row : result.covariance)
                    row[i] = nan;
                continue;
            }
            freeIndices.push_back(i);
            freeParameters.push_back(result.parameters[i]);
        }
    }

    if (freeIndices.empty())
        return;

    auto objective = [&](std::span<const double> free) {
        std::vector<double> parameters = result.parameters;
        for (std::size_t i = 0; i < freeIndices.size(); ++i)
            parameters[freeIndices[i]] = free[i];
        return naturalObjective(parameters);
    };
    auto information = math::numericalHessian(objective, freeParameters);
    auto covariance = inverse(std::move(information));
    if (covariance.empty()) {
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
        requireMethod(method, "Normal", {Method::MLE, Method::MethodOfMoments, Method::UMVU});
        Distribution distribution;
        distribution.Fit(sample, method == Method::UMVU);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::floating_point Real>
struct FitModel<StudentTRand<Real>>
{
    using Distribution = StudentTRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"degree", distribution.GetDegree(), 0.0, std::nullopt},
                {"location", distribution.GetLocation(), std::nullopt, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(parameters[0], parameters[1], parameters[2]);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Student-t", {Method::MLE});
        const double location = sampleQuantile(sample, 0.5);
        const double scale = (sampleQuantile(sample, 0.75) - sampleQuantile(sample, 0.25)) / 1.35;
        if (!(scale > 0))
            throw fit_error(FitErrorCode::WrongSample, "Student-t distribution: sample spread should be positive");
        return Distribution(5.0, location, scale);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<NoncentralTRand<Real>>
{
    using Distribution = NoncentralTRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"degree", distribution.GetDegree(), 0.0, std::nullopt},
                {"noncentrality", distribution.GetNoncentrality(), std::nullopt, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Noncentral-t", {Method::MLE});
        const auto fixedDegree = fixed.find("degree");
        const double degree = fixedDegree == fixed.end() ? 10.0 : fixedDegree->second;
        const double median = sampleQuantile(sample, 0.5);
        if (degree <= 1)
            return Distribution(degree, median);
        const long double coefficient =
            std::exp(0.5L * std::log(0.5L * degree) + std::lgamma(0.5L * (degree - 1)) - std::lgamma(0.5L * degree));
        return Distribution(degree, static_cast<double>(median / coefficient));
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<CauchyRand<Real>>
{
    using Distribution = CauchyRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"location", distribution.GetLocation(), std::nullopt, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Cauchy", {Method::MLE});
        const double scale = 0.5 * (sampleQuantile(sample, 0.75) - sampleQuantile(sample, 0.25));
        if (!(scale > 0))
            throw fit_error(FitErrorCode::WrongSample, "Cauchy distribution: sample spread should be positive");
        return Distribution(sampleQuantile(sample, 0.5), scale);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<StableRand<Real>>
{
    using Distribution = StableRand<Real>;
    using Parameterization = typename Distribution::Parameterization;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"exponent", distribution.GetExponent(), 0.1, 2.0, true, true},
                {"skewness", distribution.GetSkewness(), -1.0, 1.0, true, true},
                {"scale", distribution.GetScale(), 0.0, std::nullopt},
                {"location", distribution.GetLocation(), std::nullopt, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(parameters[0], parameters[1], parameters[2], parameters[3], Parameterization::S0);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Stable", {Method::MLE, Method::Quantile});
        const double p05 = sampleQuantile(sample, 0.05);
        const double p25 = sampleQuantile(sample, 0.25);
        const double p50 = sampleQuantile(sample, 0.50);
        const double p75 = sampleQuantile(sample, 0.75);
        const double p95 = sampleQuantile(sample, 0.95);
        const double interquartileRange = p75 - p25;
        const double outerRange = p95 - p05;
        if (!(interquartileRange > 0) || !(outerRange > 0))
            throw fit_error(FitErrorCode::WrongSample, "Stable distribution: sample spread should be positive");

        static constexpr std::array nuAlphaGrid = {2.439, 2.5, 2.6, 2.7, 2.8,  3.0,  3.2, 3.5,
                                                   4.0,   5.0, 6.0, 8.0, 10.0, 15.0, 25.0};
        static constexpr std::array nuBetaGrid = {0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0};
        static constexpr std::array alphaTable = {
            std::array{2.000, 2.000, 2.000, 2.000, 2.000, 2.000, 2.000},
            std::array{1.916, 1.924, 1.924, 1.924, 1.924, 1.924, 1.924},
            std::array{1.808, 1.813, 1.829, 1.829, 1.829, 1.829, 1.829},
            std::array{1.729, 1.730, 1.737, 1.745, 1.745, 1.745, 1.745},
            std::array{1.664, 1.663, 1.663, 1.668, 1.676, 1.676, 1.676},
            std::array{1.563, 1.560, 1.553, 1.548, 1.547, 1.547, 1.547},
            std::array{1.484, 1.480, 1.471, 1.460, 1.448, 1.438, 1.438},
            std::array{1.391, 1.386, 1.378, 1.364, 1.337, 1.318, 1.318},
            std::array{1.279, 1.273, 1.266, 1.250, 1.210, 1.184, 1.150},
            std::array{1.128, 1.121, 1.114, 1.101, 1.067, 1.027, 0.973},
            std::array{1.029, 1.021, 1.014, 1.004, 0.974, 0.935, 0.874},
            std::array{0.896, 0.892, 0.884, 0.883, 0.855, 0.823, 0.769},
            std::array{0.818, 0.812, 0.806, 0.801, 0.780, 0.756, 0.691},
            std::array{0.698, 0.695, 0.692, 0.689, 0.676, 0.656, 0.597},
            std::array{0.593, 0.590, 0.588, 0.586, 0.579, 0.563, 0.513},
        };
        static constexpr std::array betaTable = {
            std::array{0.000, 2.160, 1.000, 1.000, 1.000, 1.000, 1.000},
            std::array{0.000, 1.592, 3.390, 1.000, 1.000, 1.000, 1.000},
            std::array{0.000, 0.759, 1.800, 1.000, 1.000, 1.000, 1.000},
            std::array{0.000, 0.482, 1.048, 1.694, 1.000, 1.000, 1.000},
            std::array{0.000, 0.360, 0.760, 1.232, 2.229, 1.000, 1.000},
            std::array{0.000, 0.253, 0.518, 0.823, 1.575, 1.000, 1.000},
            std::array{0.000, 0.203, 0.410, 0.632, 1.244, 1.906, 1.000},
            std::array{0.000, 0.165, 0.332, 0.499, 0.943, 1.560, 1.000},
            std::array{0.000, 0.136, 0.271, 0.404, 0.689, 1.230, 2.195},
            std::array{0.000, 0.109, 0.216, 0.323, 0.539, 0.827, 1.917},
            std::array{0.000, 0.096, 0.190, 0.284, 0.472, 0.693, 1.759},
            std::array{0.000, 0.082, 0.163, 0.243, 0.412, 0.601, 1.596},
            std::array{0.000, 0.074, 0.147, 0.220, 0.377, 0.546, 1.482},
            std::array{0.000, 0.064, 0.128, 0.191, 0.330, 0.478, 1.362},
            std::array{0.000, 0.056, 0.112, 0.167, 0.285, 0.428, 1.274},
        };

        const double nuAlpha = outerRange / interquartileRange;
        const double nuBeta = (p95 + p05 - 2 * p50) / outerRange;
        double exponent = 2.0;
        double skewness = 0.0;
        if (nuAlpha >= nuAlphaGrid.front()) {
            exponent = bilinearInterpolation(nuAlphaGrid, nuBetaGrid, alphaTable, nuAlpha, std::fabs(nuBeta));
            const double magnitude =
                bilinearInterpolation(nuAlphaGrid, nuBetaGrid, betaTable, nuAlpha, std::fabs(nuBeta));
            skewness = std::copysign(std::clamp(magnitude, 0.0, 1.0), nuBeta);
        }
        exponent = std::clamp(exponent, 0.1, 2.0);

        const Distribution standard(exponent, skewness, 1, 0, Parameterization::S0);
        const double theoretical25 = standard.Quantile(0.25);
        const double theoretical50 = standard.Quantile(0.50);
        const double theoretical75 = standard.Quantile(0.75);
        const double scale = interquartileRange / (theoretical75 - theoretical25);
        const double location = p50 - scale * theoretical50;
        if (!(scale > 0) || !std::isfinite(location))
            throw fit_error(FitErrorCode::WrongReturn, "Stable distribution: quantile fit did not converge");

        if (method == Method::MLE) {
            exponent = std::clamp(exponent, 0.100001, 1.999999);
            skewness = std::clamp(skewness, -0.999999, 0.999999);
        }
        return Distribution(exponent, skewness, scale, location, Parameterization::S0);
    }

    static bool hasClosedForm(Method method) { return method == Method::Quantile; }
};

template <std::floating_point Real>
struct FitModel<LogStableRand<Real>>
{
    using Distribution = LogStableRand<Real>;
    using Stable = StableRand<Real>;
    using Parameterization = typename Distribution::Parameterization;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"exponent", distribution.GetExponent(), 0.1, 2.0, true, true},
                {"skewness", distribution.GetSkewness(), -1.0, 1.0, true, true},
                {"scale", distribution.GetScale(), 0.0, std::nullopt},
                {"location", distribution.GetLocation(), std::nullopt, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(parameters[0], parameters[1], parameters[2], parameters[3], Parameterization::S0);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &fixed)
    {
        if (std::ranges::any_of(sample, [](Real value) { return !(value > 0); }))
            throw fit_error(FitErrorCode::WrongSample, "Log-stable distribution: sample values should be positive");
        std::vector<Real> logSample;
        logSample.reserve(sample.size());
        for (Real value : sample)
            logSample.push_back(std::log(value));
        const Stable estimate = FitModel<Stable>::initial(logSample, method, fixed);
        return Distribution(estimate.GetExponent(), estimate.GetSkewness(), estimate.GetScale(), estimate.GetLocation(),
                            Parameterization::S0);
    }

    static bool hasClosedForm(Method method) { return FitModel<Stable>::hasClosedForm(method); }
};

template <std::floating_point Real>
struct FitModel<StableCountRand<Real>>
{
    using Distribution = StableCountRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"exponent", distribution.GetExponent(), 0.0, 1.0},
                {"scale", distribution.GetScale(), 0.0, std::nullopt},
                {"location", distribution.GetLocation(), std::nullopt, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(parameters[0], parameters[1], parameters[2]);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Stable Count", {Method::MLE});
        const auto [minimum, maximum] = std::ranges::minmax(sample);
        const double spread = static_cast<double>(maximum - minimum);
        if (!(spread > 0))
            throw fit_error(FitErrorCode::WrongSample, "Stable Count distribution: sample spread should be positive");
        const double exponent = fixed.contains("exponent") ? fixed.at("exponent") : 0.5;
        const double location =
            fixed.contains("location") ? fixed.at("location") : static_cast<double>(minimum) - 0.05 * spread;
        const double mean = static_cast<double>(sampleMeanAndVariance(sample).first);
        const Distribution standard(exponent);
        const double scale = fixed.contains("scale")
                                 ? fixed.at("scale")
                                 : std::max(static_cast<double>((mean - location) / standard.Mean()), spread * 1e-3);
        return Distribution(exponent, scale, location);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<GumbelRand<Real>>
{
    using Distribution = GumbelRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"location", distribution.GetLocation(), std::nullopt, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Gumbel", {Method::MLE});
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        const double scale = std::sqrt(6 * variance) / std::numbers::pi;
        if (!(scale > 0))
            throw fit_error(FitErrorCode::WrongSample, "Gumbel distribution: sample variance should be positive");
        return Distribution(mean - static_cast<double>(randlib::consts::egamma) * scale, scale);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<FrechetRand<Real>>
{
    using Distribution = FrechetRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), 0.0, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt},
                {"location", distribution.GetLocation(), std::nullopt, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(parameters[0], parameters[1], parameters[2]);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Frechet", {Method::MLE});
        const auto [minimum, maximum] = std::ranges::minmax(sample);
        const double spread = static_cast<double>(maximum - minimum);
        if (!(spread > 0))
            throw fit_error(FitErrorCode::WrongSample, "Frechet distribution: sample spread should be positive");
        return Distribution(1.0, spread, static_cast<double>(minimum) - 0.1 * spread);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<GeneralizedParetoRand<Real>>
{
    using Distribution = GeneralizedParetoRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), -0.5, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(parameters[0], parameters[1], 0);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Generalized Pareto", {Method::MLE});
        if (std::ranges::any_of(sample, [](Real value) { return value < 0; }))
            throw fit_error(FitErrorCode::WrongSample,
                            "Generalized Pareto distribution: exceedances should be non-negative");
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        if (!(mean > 0) || !(variance > 0))
            throw fit_error(FitErrorCode::WrongSample,
                            "Generalized Pareto distribution: sample spread should be positive");
        const double shape = std::clamp(0.5 * (1 - mean * mean / variance), -0.4, 1.0);
        return Distribution(shape, mean * (1 - shape), 0);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<GeneralizedExtremeValueRand<Real>>
{
    using Distribution = GeneralizedExtremeValueRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), std::nullopt, std::nullopt},
                {"scale", distribution.GetScale(), 0.0, std::nullopt},
                {"location", distribution.GetLocation(), std::nullopt, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(parameters[0], parameters[1], parameters[2]);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Generalized extreme-value", {Method::MLE});
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        if (!(variance > 0))
            throw fit_error(FitErrorCode::WrongSample,
                            "Generalized extreme-value distribution: sample spread should be positive");
        const double scale = std::sqrt(6 * variance) / std::numbers::pi;
        return Distribution(0, scale, mean - static_cast<double>(randlib::consts::egamma) * scale);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<FisherFRand<Real>>
{
    using Distribution = FisherFRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"degree1", static_cast<double>(distribution.GetFirstDegree()), 1.0, std::nullopt, true, false, true},
                {"degree2", static_cast<double>(distribution.GetSecondDegree()), 1.0, std::nullopt, true, false, true}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(static_cast<int>(parameters[0]), static_cast<int>(parameters[1]));
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Fisher-F", {Method::MLE});
        if (std::ranges::any_of(sample, [](Real value) { return !(value > 0); }))
            throw fit_error(FitErrorCode::WrongSample, "Fisher-F distribution: sample values should be positive");
        const auto fixedDegree1 = fixedPositiveInteger(fixed, "degree1", "Fisher-F");
        const auto fixedDegree2 = fixedPositiveInteger(fixed, "degree2", "Fisher-F");
        const int minimum1 = fixedDegree1.value_or(1);
        const int maximum1 = fixedDegree1.value_or(80);
        const int minimum2 = fixedDegree2.value_or(1);
        const int maximum2 = fixedDegree2.value_or(80);
        double bestLogLikelihood = -std::numeric_limits<double>::infinity();
        Distribution best;
        for (int degree1 = minimum1; degree1 <= maximum1; ++degree1) {
            for (int degree2 = minimum2; degree2 <= maximum2; ++degree2) {
                Distribution candidate(degree1, degree2);
                const double logLikelihood = candidate.LogLikelihoodFunction(sample);
                if (logLikelihood > bestLogLikelihood) {
                    bestLogLikelihood = logLikelihood;
                    best = std::move(candidate);
                }
            }
        }
        return best;
    }

    static bool hasClosedForm(Method) { return true; }
    static bool requiresOptimization(Method, const FixedParameterMap &) { return false; }
};

template <std::floating_point Real>
struct FitModel<InverseGammaRand<Real>>
{
    using Distribution = InverseGammaRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), 0.0, std::nullopt},
                {"rate", distribution.GetRate(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Inverse-gamma", {Method::MLE});
        if (std::ranges::any_of(sample, [](Real value) { return !(value > 0); }))
            throw fit_error(FitErrorCode::WrongSample, "Inverse-gamma distribution: sample values should be positive");
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        if (!(variance > 0))
            throw fit_error(FitErrorCode::WrongSample,
                            "Inverse-gamma distribution: sample variance should be positive");
        const double shape = 2 + mean * mean / variance;
        return Distribution(shape, mean * (shape - 1));
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<InverseGaussianRand<Real>>
{
    using Distribution = InverseGaussianRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"mean", distribution.GetMean(), 0.0, std::nullopt},
                {"shape", distribution.GetShape(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Inverse-Gaussian", {Method::MLE});
        if (std::ranges::any_of(sample, [](Real value) { return !(value > 0); }))
            throw fit_error(FitErrorCode::WrongSample,
                            "Inverse-Gaussian distribution: sample values should be positive");
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        if (!(variance > 0))
            throw fit_error(FitErrorCode::WrongSample,
                            "Inverse-Gaussian distribution: sample variance should be positive");
        return Distribution(mean, mean * mean * mean / variance);
    }

    static bool hasClosedForm(Method) { return false; }
};

template <std::floating_point Real>
struct FitModel<NoncentralChiSquaredRand<Real>>
{
    using Distribution = NoncentralChiSquaredRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"degree", distribution.GetDegree(), 0.0, std::nullopt},
                {"noncentrality", distribution.GetNoncentrality(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0], parameters[1]); }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Noncentral chi-squared", {Method::MLE});
        if (std::ranges::any_of(sample, [](Real value) { return value < 0; }))
            throw fit_error(FitErrorCode::WrongSample,
                            "Noncentral chi-squared distribution: sample values should be non-negative");
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        const double noncentrality = std::max(0.1, 0.5 * variance - mean);
        const double degree = std::max(0.1, mean - noncentrality);
        return Distribution(degree, noncentrality);
    }

    static bool hasClosedForm(Method) { return false; }
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
struct FitModel<ChiSquaredRand<Real>>
{
    using Distribution = ChiSquaredRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"degree", static_cast<double>(distribution.GetDegree()), 1.0, std::nullopt, true, false, true}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(static_cast<std::size_t>(parameters[0]));
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Chi-squared", {Method::MLE});
        if (const auto degree = fixedPositiveInteger(fixed, "degree", "chi-squared"))
            return Distribution(static_cast<std::size_t>(*degree));
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        const double searchLimit = std::clamp(std::ceil(std::max(mean, 0.5 * variance) * 4 + 20), 40.0, 1000.0);
        const int maximum = static_cast<int>(searchLimit);
        double bestLogLikelihood = -std::numeric_limits<double>::infinity();
        Distribution best;
        for (int degree = 1; degree <= maximum; ++degree) {
            Distribution candidate(static_cast<std::size_t>(degree));
            const double logLikelihood = candidate.LogLikelihoodFunction(sample);
            if (logLikelihood > bestLogLikelihood) {
                bestLogLikelihood = logLikelihood;
                best = std::move(candidate);
            }
        }
        return best;
    }

    static bool hasClosedForm(Method) { return true; }
    static bool requiresOptimization(Method, const FixedParameterMap &) { return false; }
};

template <std::floating_point Real>
struct FitModel<ErlangRand<Real>>
{
    using Distribution = ErlangRand<Real>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"shape", distribution.GetShape(), 1.0, std::nullopt, true, false, true},
                {"rate", distribution.GetRate(), 0.0, std::nullopt}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(static_cast<int>(parameters[0]), parameters[1]);
    }

    static Distribution initial(std::span<const Real> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Erlang", {Method::MLE});
        if (std::ranges::any_of(sample, [](Real value) { return !(value > 0); }))
            throw fit_error(FitErrorCode::WrongSample, "Erlang distribution: sample values should be positive");
        const auto [mean, variance] = sampleMeanAndVariance(sample);
        const auto fixedShape = fixedPositiveInteger(fixed, "shape", "Erlang");
        const auto fixedRate = fixed.find("rate");
        if (fixedRate != fixed.end() && !(fixedRate->second > 0))
            throw std::invalid_argument("Fixed Erlang rate should be positive");
        const double estimate = variance > 0 ? std::clamp(std::round(mean * mean / variance), 1.0, 245.0) : 1.0;
        const int maximum = fixedShape ? *fixedShape : std::clamp(static_cast<int>(4 * estimate + 20), 40, 1000);
        const int minimum = fixedShape ? *fixedShape : 1;
        double bestLogLikelihood = -std::numeric_limits<double>::infinity();
        Distribution best;
        for (int shape = minimum; shape <= maximum; ++shape) {
            const double rate = fixedRate == fixed.end() ? shape / mean : fixedRate->second;
            Distribution candidate(shape, rate);
            const double logLikelihood = candidate.LogLikelihoodFunction(sample);
            if (logLikelihood > bestLogLikelihood) {
                bestLogLikelihood = logLikelihood;
                best = std::move(candidate);
            }
        }
        return best;
    }

    static bool hasClosedForm(Method) { return true; }
    static bool requiresOptimization(Method, const FixedParameterMap &) { return false; }
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
        const auto [minimum, maximum] = std::ranges::minmax(sample);
        const double spread = static_cast<double>(maximum - minimum);
        if (!(spread > 0))
            throw fit_error(FitErrorCode::WrongSample, "Levy distribution: sample spread should be positive");
        const auto fixedLocation = fixed.find("location");
        const double location =
            fixedLocation == fixed.end() ? static_cast<double>(minimum) - 0.1 * spread : fixedLocation->second;
        Distribution distribution(location, spread);
        distribution.FitScale(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return false; }
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
        requireMethod(method, "Poisson", {Method::MLE, Method::MethodOfMoments, Method::UMVU});
        Distribution distribution;
        distribution.Fit(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
};

template <>
struct FitModel<BernoulliRand>
{
    using Distribution = BernoulliRand;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"probability", distribution.GetProbability(), 0.0, 1.0, true, true}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0]); }

    static Distribution initial(std::span<const int> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Bernoulli", {Method::MLE, Method::MethodOfMoments});
        if (std::ranges::any_of(sample, [](int value) { return value != 0 && value != 1; }))
            throw fit_error(FitErrorCode::WrongSample, "Bernoulli distribution: sample values should be zero or one");
        const double sum = std::accumulate(sample.begin(), sample.end(), 0.0);
        return Distribution(sum / sample.size());
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::integral Int>
struct FitModel<GeometricRand<Int>>
{
    using Distribution = GeometricRand<Int>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"probability", distribution.GetProbability(), 0.0, 1.0, false, true}};
    }

    static Distribution make(std::span<const double> parameters) { return Distribution(parameters[0]); }

    static Distribution initial(std::span<const Int> sample, Method method, const FixedParameterMap &)
    {
        requireMethod(method, "Geometric", {Method::MLE, Method::MethodOfMoments});
        if (std::ranges::any_of(sample, [](Int value) { return value < 0; }))
            throw fit_error(FitErrorCode::WrongSample, "Geometric distribution: sample values should be non-negative");
        const double sum = std::accumulate(sample.begin(), sample.end(), 0.0);
        return Distribution(1 / (1 + sum / sample.size()));
    }

    static bool hasClosedForm(Method) { return true; }
};

template <std::integral Int>
struct FitModel<BinomialRand<Int>>
{
    using Distribution = BinomialRand<Int>;

    static std::vector<FitParameter> parameters(const Distribution &distribution)
    {
        return {{"number", static_cast<double>(distribution.GetNumber()), 1.0,
                 static_cast<double>(std::numeric_limits<int>::max()), true, true, true},
                {"probability", distribution.GetProbability(), 0.0, 1.0, true, true}};
    }

    static Distribution make(std::span<const double> parameters)
    {
        return Distribution(static_cast<int>(parameters[0]), parameters[1]);
    }

    static Distribution initial(std::span<const Int> sample, Method method, const FixedParameterMap &fixed)
    {
        requireMethod(method, "Binomial", {Method::MLE});
        const auto number = fixedPositiveInteger(fixed, "number", "binomial");
        if (!number)
            throw fit_error(FitErrorCode::NotApplicable,
                            "Binomial distribution: fitting currently requires a fixed number");
        Distribution distribution(*number);
        distribution.FitProbability(sample);
        return distribution;
    }

    static bool hasClosedForm(Method) { return true; }
    static bool requiresOptimization(Method, const FixedParameterMap &) { return false; }
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
        return {{"minimum", static_cast<double>(distribution.MinValue()),
                 static_cast<double>(std::numeric_limits<Int>::lowest()), static_cast<double>(distribution.MaxValue()),
                 true, false, true},
                {"maximum", static_cast<double>(distribution.MaxValue()), static_cast<double>(distribution.MinValue()),
                 static_cast<double>(std::numeric_limits<Int>::max()), false, true, true}};
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

template <class Distribution, class LogLikelihood>
FitResult<Distribution> fitLikelihood(std::span<const typename Distribution::value_type> initializationSample,
                                      FitOptions options, const LogLikelihood &logLikelihood,
                                      bool forceOptimization = false)
{
    if (initializationSample.empty())
        throw fit_error(FitErrorCode::TooFewElements, "Fitting requires a non-empty sample");
    if constexpr (std::floating_point<typename Distribution::value_type>) {
        if (std::ranges::any_of(initializationSample, [](auto value) { return !std::isfinite(value); }))
            throw fit_error(FitErrorCode::WrongSample, "Fitting requires finite sample values");
    }
    if (!(options.tolerance > 0) || options.max_iterations <= 0)
        throw std::invalid_argument("Fitting tolerance and iteration limit should be positive");
    if (options.method != Method::MLE && !options.fixed.empty())
        throw fit_error(FitErrorCode::NotApplicable, "Fixed parameters are currently supported only for MLE");

    Distribution initial = FitModel<Distribution>::initial(initializationSample, options.method, options.fixed);
    auto specifications = FitModel<Distribution>::parameters(initial);
    std::vector<double> parameters;
    parameters.reserve(specifications.size());
    std::vector<bool> fixed(specifications.size());
    for (std::size_t i = 0; i < specifications.size(); ++i) {
        auto found = options.fixed.find(specifications[i].name);
        if (found != options.fixed.end()) {
            if (!withinBounds(found->second, specifications[i]))
                throw std::invalid_argument("Fixed parameter is outside its bounds: " +
                                            std::string(specifications[i].name));
            specifications[i].value = found->second;
            fixed[i] = true;
        }
        parameters.push_back(specifications[i].value);
    }
    for (const auto &[name, value] : options.fixed) {
        (void)value;
        const std::string &fixedName = name;
        if (std::ranges::none_of(specifications,
                                 [&fixedName](const FitParameter &parameter) { return parameter.name == fixedName; }))
            throw std::invalid_argument("Unknown fixed parameter: " + name);
    }

    const bool optimize = forceOptimization || requiresOptimization<Distribution>(options.method, options.fixed);
    std::vector<std::size_t> freeIndices;
    std::vector<double> unconstrained;
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        if (!fixed[i]) {
            if (specifications[i].integral && optimize)
                throw fit_error(FitErrorCode::NotApplicable,
                                "Free integer parameters are not supported by continuous optimization");
            freeIndices.push_back(i);
            unconstrained.push_back(toUnconstrained(parameters[i], specifications[i]));
        }
    }

    auto naturalObjective = [&](std::span<const double> candidate) {
        return negativeLogLikelihood<Distribution>(candidate, logLikelihood);
    };
    math::MinimizeResult<double> optimization;
    optimization.x = unconstrained;
    optimization.converged = true;
    if (!freeIndices.empty() && optimize) {
        auto objective = [&](std::span<const double> free) {
            std::vector<double> candidate = parameters;
            for (std::size_t i = 0; i < free.size(); ++i)
                candidate[freeIndices[i]] = fromUnconstrained(free[i], specifications[freeIndices[i]]);
            return naturalObjective(candidate);
        };
        if (forceOptimization)
            optimization = math::nelderMead(objective, unconstrained, options.tolerance, options.max_iterations);
        else {
            optimization = math::bfgs(objective, unconstrained, options.tolerance, options.max_iterations);
            if (!optimization.converged)
                optimization =
                    math::nelderMead(objective, std::move(optimization.x), options.tolerance, options.max_iterations);
        }
        for (std::size_t i = 0; i < optimization.x.size(); ++i)
            parameters[freeIndices[i]] = fromUnconstrained(optimization.x[i], specifications[freeIndices[i]]);
    }

    Distribution distribution = FitModel<Distribution>::make(parameters);
    const double fittedLogLikelihood = logLikelihood(distribution);
    const std::size_t fittedParameters = freeIndices.size();
    FitResult<Distribution> result(std::move(distribution));
    result.log_likelihood = fittedLogLikelihood;
    result.aic = 2 * fittedParameters - 2 * fittedLogLikelihood;
    result.bic =
        std::log(static_cast<double>(initializationSample.size())) * fittedParameters - 2 * fittedLogLikelihood;
    result.parameters = std::move(parameters);
    result.fixed = std::move(fixed);
    result.converged = optimization.converged;
    result.iterations = optimization.iterations;
    result.method = options.method;
    for (const auto &parameter : specifications)
        result.parameter_names.emplace_back(parameter.name);
    addUncertainty(result, specifications, naturalObjective);
    return result;
}

} // namespace detail

template <class Distribution>
    requires detail::FittableModel<Distribution>
FitResult<Distribution> fit(std::span<const typename Distribution::value_type> sample, FitOptions options = {})
{
    const auto logLikelihood = [sample](const Distribution &distribution) {
        return distribution.LogLikelihoodFunction(sample);
    };
    return detail::fitLikelihood<Distribution>(sample, std::move(options), logLikelihood);
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

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif

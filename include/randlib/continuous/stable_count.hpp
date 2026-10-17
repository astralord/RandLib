#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>

#include <randlib/continuous/stable.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/math/constants.hpp>
#include <randlib/math/quadrature.hpp>

namespace randlib
{

/**
 * @brief Stable Count distribution
 *
 * The standard density is
 *
 * \f[
 * N_\alpha(x) = \frac{L_\alpha(1/x)}{\Gamma(1/\alpha+1)x}, \qquad x>0,
 * \f]
 *
 * where the positive stable density has Laplace transform
 * \f$\mathrm{E}[\exp(-sS)]=\exp(-s^\alpha)\f$. In RandLib's S1 convention this is
 * `StableRand(alpha, 1, pow(cos(pi * alpha / 2), 1 / alpha), 0, S1)`.
 * `StableCountRand(alpha, scale, location)` is the law of `location + scale * X`
 * for standard Stable Count `X`.
 *
 * At \f$\alpha=1/2\f$, the standard law is exactly Gamma(shape=3/2, scale=4).
 * As \f$\alpha\to1\f$, it converges weakly to a point mass at one; the endpoint
 * \f$\alpha=1\f$ itself is deliberately excluded from this continuous distribution.
 */
template <typename RealType>
class StableCountRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    double alpha = 0.5;
    double theta = 1;
    double nu = 0;
    double logTheta = 0;
    double logNormalization = std::lgamma(3.0);
    StableRand<RealType> positiveStable{0.5, 1.0, 0.5, 0.0};

    [[nodiscard]] double standardLogDensity(double x) const;
    [[nodiscard]] double wrightSeries(double x, bool integrated) const;
    [[nodiscard]] double directStandardCdf(double x) const;
    [[nodiscard]] double directStandardSurvival(double x) const;
    [[nodiscard]] long double standardMoment(long double order) const;

    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

public:
    StableCountRand(double exponent = 0.5, double scale = 1, double location = 0)
    {
        SetParameters(exponent, scale, location);
    }

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return static_cast<RealType>(nu); }
    [[nodiscard]] RealType MaxValue() const override { return std::numeric_limits<RealType>::infinity(); }

    void SetParameters(double exponent, double scale = 1, double location = 0);
    void SetExponent(double exponent) { SetParameters(exponent, theta, nu); }
    void SetScale(double scale) { SetParameters(alpha, scale, nu); }
    void SetLocation(double location) { SetParameters(alpha, theta, location); }

    [[nodiscard]] double GetExponent() const noexcept { return alpha; }
    [[nodiscard]] double GetScale() const noexcept { return theta; }
    [[nodiscard]] double GetLocation() const noexcept { return nu; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    /** @return the standard raw moment E[X^order], defined for order > -1 */
    [[nodiscard]] long double StandardMoment(long double order) const { return standardMoment(order); }
    /** @return the nonstandard raw integer moment E[X^order] */
    [[nodiscard]] long double Moment(int order) const;
    [[nodiscard]] long double Mean() const override { return Moment(1); }
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] long double SecondMoment() const override { return Moment(2); }
    [[nodiscard]] long double ThirdMoment() const override { return Moment(3); }
    [[nodiscard]] long double FourthMoment() const override { return Moment(4); }
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
};

template <typename RealType>
void StableCountRand<RealType>::SetParameters(double exponent, double scale, double location)
{
    if (!(exponent > 0 && exponent < 1) || !std::isfinite(exponent))
        throw std::invalid_argument("Stable Count distribution: exponent should be in the interval (0, 1)");
    if (!(scale > 0) || !std::isfinite(scale))
        throw std::invalid_argument("Stable Count distribution: scale should be finite and positive");
    if (!std::isfinite(location))
        throw std::invalid_argument("Stable Count distribution: location should be finite");

    alpha = exponent;
    theta = scale;
    nu = location;
    logTheta = std::log(theta);
    logNormalization = std::lgamma(1 / alpha + 1);
    const double stableScale = std::exp(std::log(std::cos(randlib::consts::half_pi * alpha)) / alpha);
    positiveStable = StableRand<RealType>(alpha, 1.0, stableScale, 0.0, StableRand<RealType>::Parameterization::S1);
}

template <typename RealType>
std::string StableCountRand<RealType>::Name() const
{
    return "Stable Count(" + this->toStringWithPrecision(alpha) + ", " + this->toStringWithPrecision(theta) + ", " +
           this->toStringWithPrecision(nu) + ")";
}

template <typename RealType>
double StableCountRand<RealType>::standardLogDensity(double x) const
{
    if (!(x > 0) || !std::isfinite(x))
        return -INFINITY;
    if (alpha == 0.5)
        return 0.5 * std::log(x) - 0.25 * x - std::log(4.0) - 0.5 * static_cast<double>(randlib::consts::lnpi);
    const double series = wrightSeries(x, false);
    if (series > 0 && std::isfinite(series))
        return std::log(series);
    return positiveStable.logf(static_cast<RealType>(1 / x)) - logNormalization - std::log(x);
}

template <typename RealType>
double StableCountRand<RealType>::wrightSeries(double x, bool integrated) const
{
    if (!(x > 0) || !std::isfinite(x) || alpha == 0.5)
        return NAN;
    const long double exponent = alpha;
    const long double logX = std::log(static_cast<long double>(x));
    const long double logPi = std::log(randlib::consts::pi);
    long double sum = 0;
    long double correction = 0;
    long double largest = 0;
    int largestIndex = 0;
    int smallTerms = 0;
    for (int k = 1; k <= 4096; ++k) {
        const long double ak = exponent * k;
        const long double sine = std::sin(randlib::consts::pi * ak);
        if (sine == 0)
            continue;
        long double logMagnitude = std::lgamma(1 + ak) - std::lgamma(static_cast<long double>(k) + 1) + ak * logX +
                                   std::log(std::fabs(sine)) - logPi - static_cast<long double>(logNormalization);
        if (integrated)
            logMagnitude += logX - std::log(ak + 1);
        if (logMagnitude > 15)
            return NAN;
        const long double magnitude = std::exp(logMagnitude);
        const long double sign = ((k % 2 == 0) ? -1 : 1) * std::copysign(1.0L, sine);
        const long double term = sign * magnitude;
        const long double adjusted = term - correction;
        const long double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
        if (magnitude > largest) {
            largest = magnitude;
            largestIndex = k;
            smallTerms = 0;
        } else if (k > largestIndex + 8 &&
                   magnitude <= 8 * std::numeric_limits<long double>::epsilon() * std::max(std::fabs(sum), 1.0L)) {
            ++smallTerms;
            if (smallTerms >= 24)
                return (sum >= 0) ? static_cast<double>(sum) : NAN;
        } else {
            smallTerms = 0;
        }
    }
    return NAN;
}

template <typename RealType>
double StableCountRand<RealType>::directStandardCdf(double x) const
{
    if (!(x > 0))
        return 0;
    if (std::isinf(x))
        return 1;
    const double series = wrightSeries(x, true);
    if (series >= 0 && series <= 1 && std::isfinite(series))
        return series;
    const double normalization = std::exp(logNormalization);
    return std::clamp(math::integrate<double>(
                          [this, normalization](double t) {
                              return positiveStable.f(static_cast<RealType>(t)) / (normalization * t);
                          },
                          1 / x, INFINITY, 2e-13, 1024),
                      0.0, 1.0);
}

template <typename RealType>
double StableCountRand<RealType>::directStandardSurvival(double x) const
{
    if (!(x > 0))
        return 1;
    if (std::isinf(x))
        return 0;
    const double normalization = std::exp(logNormalization);
    return std::clamp(math::integrate<double>(
                          [this, normalization](double t) {
                              return positiveStable.f(static_cast<RealType>(t)) / (normalization * t);
                          },
                          0.0, 1 / x, 2e-13, 1024),
                      0.0, 1.0);
}

template <typename RealType>
double StableCountRand<RealType>::logf(const RealType &x) const
{
    const double standardized = (static_cast<double>(x) - nu) / theta;
    return standardLogDensity(standardized) - logTheta;
}

template <typename RealType>
double StableCountRand<RealType>::f(const RealType &x) const
{
    return std::exp(logf(x));
}

template <typename RealType>
double StableCountRand<RealType>::F(const RealType &x) const
{
    const double standardized = (static_cast<double>(x) - nu) / theta;
    if (!(standardized > 0))
        return 0;
    if (alpha == 0.5) {
        static const double logShape = std::log(1.5);
        static const double logGammaShape = std::lgamma(1.5);
        return math::pgamma(1.5, 0.25 * standardized, logShape, logGammaShape);
    }
    const double lower = directStandardCdf(standardized);
    return (lower <= 0.95) ? lower : 1 - directStandardSurvival(standardized);
}

template <typename RealType>
double StableCountRand<RealType>::S(const RealType &x) const
{
    const double standardized = (static_cast<double>(x) - nu) / theta;
    if (!(standardized > 0))
        return 1;
    if (alpha == 0.5) {
        static const double logShape = std::log(1.5);
        static const double logGammaShape = std::lgamma(1.5);
        return math::qgamma(1.5, 0.25 * standardized, logShape, logGammaShape);
    }
    const double lower = directStandardCdf(standardized);
    return (lower <= 0.95) ? 1 - lower : directStandardSurvival(standardized);
}

template <typename RealType>
RealType StableCountRand<RealType>::quantileImpl(double p) const
{
    const long double mean = standardMoment(1);
    const RealType guess =
        std::isfinite(mean) ? static_cast<RealType>(nu + theta * mean) : std::numeric_limits<RealType>::quiet_NaN();
    return ContinuousDistribution<RealType>::quantileImpl(p, guess);
}

template <typename RealType>
RealType StableCountRand<RealType>::quantileImpl1m(double p) const
{
    const long double mean = standardMoment(1);
    const RealType guess =
        std::isfinite(mean) ? static_cast<RealType>(nu + theta * mean) : std::numeric_limits<RealType>::quiet_NaN();
    return ContinuousDistribution<RealType>::quantileImpl1m(p, guess);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType StableCountRand<RealType>::operator()(G &g) const
{
    if (alpha == 0.5) {
        const RealType exponential = detail::standardExponential<RealType>(g);
        const RealType normal = detail::standardNormal<RealType>(g);
        const RealType standard = 4 * (exponential + 0.5 * normal * normal);
        return static_cast<RealType>(nu) + static_cast<RealType>(theta) * standard;
    }
    return this->Quantile(detail::uniform01<double>(g));
}

template <typename RealType>
void StableCountRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    for (RealType &value : outputData)
        value = (*this)(g);
}

template <typename RealType>
long double StableCountRand<RealType>::standardMoment(long double order) const
{
    if (!(order > -1) || !std::isfinite(order))
        throw std::invalid_argument("Stable Count distribution: standard moment order should be finite and above -1");
    const long double inverseAlpha = 1 / static_cast<long double>(alpha);
    const long double logarithm =
        std::lgamma((order + 1) * inverseAlpha) - std::lgamma(order + 1) - std::lgamma(inverseAlpha);
    return std::exp(logarithm);
}

template <typename RealType>
long double StableCountRand<RealType>::Moment(int order) const
{
    if (order < 0)
        throw std::invalid_argument("Stable Count distribution: raw moment order should be non-negative");
    long double sum = 0;
    long double binomial = 1;
    for (int k = 0; k <= order; ++k) {
        sum += binomial * std::pow(static_cast<long double>(nu), order - k) *
               std::pow(static_cast<long double>(theta), k) * standardMoment(k);
        binomial *= static_cast<long double>(order - k) / (k + 1);
    }
    return sum;
}

template <typename RealType>
long double StableCountRand<RealType>::Variance() const
{
    const long double first = standardMoment(1);
    const long double second = standardMoment(2);
    return static_cast<long double>(theta) * theta * (second - first * first);
}

template <typename RealType>
long double StableCountRand<RealType>::Skewness() const
{
    const long double first = standardMoment(1);
    const long double second = standardMoment(2);
    const long double third = standardMoment(3);
    const long double variance = second - first * first;
    return (third - 3 * first * second + 2 * first * first * first) / std::pow(variance, 1.5L);
}

template <typename RealType>
long double StableCountRand<RealType>::ExcessKurtosis() const
{
    const long double first = standardMoment(1);
    const long double second = standardMoment(2);
    const long double third = standardMoment(3);
    const long double fourth = standardMoment(4);
    const long double variance = second - first * first;
    const long double centralFourth =
        fourth - 4 * first * third + 6 * first * first * second - 3 * first * first * first * first;
    return centralFourth / (variance * variance) - 3;
}

} // namespace randlib

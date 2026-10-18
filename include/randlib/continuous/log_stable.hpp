#pragma once

#include <cmath>
#include <complex>

#include <randlib/continuous/stable.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/math/constants.hpp>
#include <randlib/math/optimize.hpp>

namespace randlib
{

/**
 * @brief Exponential transform of a stable distribution
 *
 * Notation: Y ~ Log-Stable(alpha, beta, gamma, mu), where log(Y) is stable.
 *
 * Positive integer moments exist for beta = -1 and for the Gaussian alpha = 2
 * case. The characteristic function is intentionally not implemented.
 */
template <typename RealType>
class LogStableRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return std::exp(stable(g)); }

    StableRand<RealType> stable{};

    [[nodiscard]] bool hasMoments() const noexcept;
    [[nodiscard]] long double logRawMoment(int order) const;

public:
    using Parameterization = typename StableRand<RealType>::Parameterization;

    LogStableRand(double exponent = 2, double skewness = 0, double scale = 1, double location = 0,
                  Parameterization parameterization = Parameterization::S1)
        : stable(exponent, skewness, scale, location, parameterization)
    {
    }

    [[nodiscard]] std::string Name() const override { return "Log-" + stable.Name(); }
    [[nodiscard]] Support SupportType() const override;
    [[nodiscard]] RealType MinValue() const override;
    [[nodiscard]] RealType MaxValue() const override;

    void SetParameters(double exponent, double skewness, double scale = 1, double location = 0,
                       Parameterization parameterization = Parameterization::S1);
    void SetExponent(double exponent) { stable.SetExponent(exponent); }
    void SetSkewness(double skewness) { stable.SetSkewness(skewness); }
    void SetScale(double scale) { stable.SetScale(scale); }
    void SetLocation(double location) { stable.SetLocation(location); }
    void SetParameterization(Parameterization parameterization) noexcept
    {
        stable.SetParameterization(parameterization);
    }

    [[nodiscard]] double GetExponent() const noexcept { return stable.GetExponent(); }
    [[nodiscard]] double GetSkewness() const noexcept { return stable.GetSkewness(); }
    [[nodiscard]] double GetScale() const noexcept { return stable.GetScale(); }
    [[nodiscard]] double GetLocation() const noexcept { return stable.GetLocation(); }
    [[nodiscard]] Parameterization GetParameterization() const noexcept { return stable.GetParameterization(); }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const
    {
        return std::exp(stable(g));
    }

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
    [[nodiscard]] RealType Median() const override { return std::exp(stable.Median()); }
    [[nodiscard]] RealType Mode() const override;

private:
    RealType quantileImpl(double p) const override { return std::exp(stable.Quantile(p)); }
    RealType quantileImpl1m(double p) const override { return std::exp(stable.Quantile1m(p)); }
    std::complex<double> CFImpl(double t) const override;
};

template <typename RealType>
void LogStableRand<RealType>::SetParameters(double exponent, double skewness, double scale, double location,
                                            Parameterization parameterization)
{
    stable = StableRand<RealType>(exponent, skewness, scale, location, parameterization);
}

template <typename RealType>
Support LogStableRand<RealType>::SupportType() const
{
    return std::isfinite(stable.MaxValue()) ? Support::Finite : Support::RightSemiFinite;
}

template <typename RealType>
RealType LogStableRand<RealType>::MinValue() const
{
    const RealType minimum = stable.MinValue();
    return std::isfinite(minimum) ? std::exp(minimum) : RealType{0};
}

template <typename RealType>
RealType LogStableRand<RealType>::MaxValue() const
{
    const RealType maximum = stable.MaxValue();
    return std::isfinite(maximum) ? std::exp(maximum) : static_cast<RealType>(INFINITY);
}

template <typename RealType>
double LogStableRand<RealType>::f(const RealType &x) const
{
    return (x >= 0) ? std::exp(logf(x)) : 0;
}

template <typename RealType>
double LogStableRand<RealType>::logf(const RealType &x) const
{
    if (x < 0)
        return -INFINITY;
    if (x == 0)
        return (stable.GetExponent() < 2 && stable.GetSkewness() < 1) ? INFINITY : -INFINITY;
    const RealType logX = std::log(x);
    return stable.logf(logX) - logX;
}

template <typename RealType>
double LogStableRand<RealType>::F(const RealType &x) const
{
    return (x > 0) ? stable.F(std::log(x)) : 0;
}

template <typename RealType>
double LogStableRand<RealType>::S(const RealType &x) const
{
    return (x > 0) ? stable.S(std::log(x)) : 1;
}

template <typename RealType>
bool LogStableRand<RealType>::hasMoments() const noexcept
{
    return stable.GetExponent() == 2 || stable.GetSkewness() == -1;
}

template <typename RealType>
long double LogStableRand<RealType>::logRawMoment(int order) const
{
    StableRand<RealType> s1 = stable;
    s1.SetParameterization(Parameterization::S1);
    const long double alpha = stable.GetExponent();
    const long double momentOrder = order;
    const long double scaleOrder = static_cast<long double>(stable.GetScale()) * momentOrder;
    const long double location = s1.GetLocation();
    if (alpha == 1)
        return momentOrder * location + randlib::consts::two_inv_pi * scaleOrder * std::log(momentOrder);
    return momentOrder * location - std::pow(scaleOrder, alpha) / std::cos(randlib::consts::half_pi * alpha);
}

template <typename RealType>
long double LogStableRand<RealType>::Mean() const
{
    return hasMoments() ? std::exp(logRawMoment(1)) : INFINITY;
}

template <typename RealType>
long double LogStableRand<RealType>::Variance() const
{
    if (!hasMoments())
        return INFINITY;
    const long double first = logRawMoment(1);
    return std::exp(2 * first) * std::expm1(logRawMoment(2) - 2 * first);
}

template <typename RealType>
long double LogStableRand<RealType>::Skewness() const
{
    if (!hasMoments())
        return NAN;
    if (stable.GetExponent() == 2) {
        const long double varianceRatio = std::expm1(2 * stable.GetScale() * stable.GetScale());
        return (varianceRatio + 3) * std::sqrt(varianceRatio);
    }
    const long double first = logRawMoment(1);
    const long double second = logRawMoment(2) - 2 * first;
    const long double third = logRawMoment(3) - 3 * first;
    const long double variance = std::expm1(second);
    return (std::expm1(third) - 3 * std::expm1(second)) / std::pow(variance, 1.5L);
}

template <typename RealType>
long double LogStableRand<RealType>::ExcessKurtosis() const
{
    if (!hasMoments())
        return NAN;
    if (stable.GetExponent() == 2) {
        const long double logVariance = 2 * stable.GetScale() * stable.GetScale();
        return std::expm1(4 * logVariance) + 2 * std::expm1(3 * logVariance) + 3 * std::expm1(2 * logVariance);
    }
    const long double first = logRawMoment(1);
    const long double second = logRawMoment(2) - 2 * first;
    const long double third = logRawMoment(3) - 3 * first;
    const long double fourth = logRawMoment(4) - 4 * first;
    const long double variance = std::expm1(second);
    return (std::expm1(fourth) - 4 * std::expm1(third) + 6 * std::expm1(second)) / (variance * variance) - 3;
}

template <typename RealType>
RealType LogStableRand<RealType>::Mode() const
{
    if (stable.GetExponent() < 2 && stable.GetSkewness() < 1)
        return 0;
    if (stable.GetExponent() == 2) {
        StableRand<RealType> s1 = stable;
        s1.SetParameterization(Parameterization::S1);
        return std::exp(static_cast<RealType>(s1.GetLocation() - 2 * stable.GetScale() * stable.GetScale()));
    }
    RealType mode = stable.Median();
    math::findMin<RealType>([this](RealType logValue) { return -stable.logf(logValue) + logValue; }, mode, mode);
    return std::exp(mode);
}

template <typename RealType>
std::complex<double> LogStableRand<RealType>::CFImpl(double t) const
{
    return this->integrateCharacteristicFunction(t);
}

} // namespace randlib

#pragma once

#include <array>
#include <functional>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/exponential_family.hpp>
#include <randlib/continuous/gamma.hpp>

namespace randlib
{

/**
 * @brief The ExponentialRand class <BR>
 * Exponential distribution
 *
 * f(x | β) = β exp(-βx)
 *
 * Notation: X ~ Exp(β)
 *
 * Related distributions: <BR>
 * X ~ Γ(1, β)
 */
template <typename RealType>
class ExponentialRand final : public FreeRateGammaDistribution<RealType>, public ExponentialFamily<RealType, double>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

public:
    explicit ExponentialRand(double rate = 1) : FreeRateGammaDistribution<RealType>(1, rate) {}

    /// @return Exp(β) with rate β
    [[nodiscard]] static ExponentialRand FromRate(double rate) { return ExponentialRand(rate); }

    /// @return Exp(1/θ) with scale (mean) θ
    [[nodiscard]] static ExponentialRand FromScale(double scale)
    {
        if (!(scale > 0.0))
            throw std::invalid_argument("Exponential distribution: scale should be positive");
        return ExponentialRand(1.0 / scale);
    }
    [[nodiscard]] std::string Name() const override;

    double SufficientStatistic(RealType x) const override;
    double SourceParameters() const override;
    double SourceToNatural(double rate) const override;
    double LogNormalizer(double theta) const override;
    double LogNormalizerGradient(double theta) const override;
    double CarrierMeasure(RealType) const override;
    double CrossEntropyAdjusted(double rate) const override;
    double EntropyAdjusted() const override;

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);

private:
    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Entropy() const;
    [[nodiscard]] long double Moment(int n) const;
    [[nodiscard]] long double ThirdMoment() const override { return Moment(3); }
    [[nodiscard]] long double FourthMoment() const override { return Moment(4); }
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template <typename RealType>
std::string ExponentialRand<RealType>::Name() const
{
    return "Exponential(" + this->toStringWithPrecision(this->GetRate()) + ")";
}

template <typename RealType>
double ExponentialRand<RealType>::SufficientStatistic(RealType x) const
{
    return x;
}

template <typename RealType>
double ExponentialRand<RealType>::SourceParameters() const
{
    return this->beta;
}

template <typename RealType>
double ExponentialRand<RealType>::SourceToNatural(double rate) const
{
    return -rate;
}

template <typename RealType>
double ExponentialRand<RealType>::LogNormalizer(double thetaP) const
{
    return -std::log(-thetaP);
}

template <typename RealType>
double ExponentialRand<RealType>::LogNormalizerGradient(double thetaP) const
{
    return -1.0 / thetaP;
}

template <typename RealType>
double ExponentialRand<RealType>::CarrierMeasure(RealType) const
{
    return 0.0;
}

template <typename RealType>
double ExponentialRand<RealType>::CrossEntropyAdjusted(double rate) const
{
    return rate * this->theta - std::log(rate);
}

template <typename RealType>
double ExponentialRand<RealType>::EntropyAdjusted() const
{
    return 1.0 - this->logBeta;
}

template <typename RealType>
double ExponentialRand<RealType>::f(const RealType &x) const
{
    return (x < 0.0) ? 0.0 : this->beta * std::exp(-this->beta * x);
}

template <typename RealType>
double ExponentialRand<RealType>::logf(const RealType &x) const
{
    return (x < 0.0) ? -INFINITY : this->logBeta - this->beta * x;
}

template <typename RealType>
double ExponentialRand<RealType>::F(const RealType &x) const
{
    return (x > 0.0) ? -std::expm1l(-this->beta * x) : 0.0;
}

template <typename RealType>
double ExponentialRand<RealType>::S(const RealType &x) const
{
    return (x > 0.0) ? std::exp(-this->beta * x) : 1.0;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType ExponentialRand<RealType>::operator()(G &g) const
{
    return this->theta * StandardVariate(g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType ExponentialRand<RealType>::StandardVariate(G &g)
{
    return detail::standardExponential<RealType>(g);
}

template <typename RealType>
std::complex<double> ExponentialRand<RealType>::CFImpl(double t) const
{
    return 1.0 / std::complex<double>(1.0, -this->theta * t);
}

template <typename RealType>
long double ExponentialRand<RealType>::Entropy() const
{
    return this->EntropyAdjusted();
}

template <typename RealType>
long double ExponentialRand<RealType>::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return std::exp(math::lfact(n) - n * this->logBeta);
}

} // namespace randlib

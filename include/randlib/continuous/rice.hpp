#pragma once

#include <span>

#include <randlib/core/continuous.hpp>
#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>

namespace randlib
{

/**
 * @brief Rice distribution
 *
 * Notation: X ~ Rice(ν, σ)
 *
 * If Y and Z are independent standard normal variables, then
 * X = hypot(ν + σY, σZ).
 */
template <typename RealType>
class RiceRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    double nu = 0;
    double sigma = 1;
    double halfNoncentrality = 0;
    double logSigma = 0;

public:
    explicit RiceRand(double offset = 0, double scale = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double offset, double scale);
    [[nodiscard]] double GetOffset() const noexcept { return nu; }
    [[nodiscard]] double GetScale() const noexcept { return sigma; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

namespace randlib
{

template <typename RealType>
RiceRand<RealType>::RiceRand(double offset, double scale)
{
    SetParameters(offset, scale);
}

template <typename RealType>
std::string RiceRand<RealType>::Name() const
{
    return "Rice(" + this->toStringWithPrecision(nu) + ", " + this->toStringWithPrecision(sigma) + ")";
}

template <typename RealType>
void RiceRand<RealType>::SetParameters(double offset, double scale)
{
    if (!(offset >= 0) || !std::isfinite(offset))
        throw std::invalid_argument("Rice distribution: offset should be finite and non-negative");
    if (!(scale > 0) || !std::isfinite(scale))
        throw std::invalid_argument("Rice distribution: scale should be finite and positive");
    nu = offset;
    sigma = scale;
    const double standardizedOffset = nu / sigma;
    halfNoncentrality = 0.5 * standardizedOffset * standardizedOffset;
    logSigma = std::log(sigma);
}

template <typename RealType>
double RiceRand<RealType>::f(const RealType &x) const
{
    return (x > 0 && std::isfinite(x)) ? std::exp(logf(x)) : 0;
}

template <typename RealType>
double RiceRand<RealType>::logf(const RealType &x) const
{
    if (!(x > 0) || !std::isfinite(x))
        return -INFINITY;
    const double standardizedDifference = (static_cast<double>(x) - nu) / sigma;
    const double z = (static_cast<double>(x) / sigma) * (nu / sigma);
    double scaledLogBessel;
    if (z > 50) {
        const double inverseZ = 1 / z;
        const double series =
            1 + inverseZ * (1.0 / 8 + inverseZ * (9.0 / 128 + inverseZ * (225.0 / 3072 + inverseZ * 11025.0 / 98304)));
        scaledLogBessel = -0.5 * std::log(2 * randlib::consts::pi * z) + std::log(series);
    } else {
        scaledLogBessel = math::logBesselI(0.0, z) - z;
    }
    return std::log(static_cast<double>(x)) - 2 * logSigma - 0.5 * standardizedDifference * standardizedDifference +
           scaledLogBessel;
}

template <typename RealType>
double RiceRand<RealType>::F(const RealType &x) const
{
    if (x <= 0)
        return 0;
    if (std::isinf(x))
        return 1;
    const double standardizedX = static_cast<double>(x) / sigma;
    const double halfX = 0.5 * standardizedX * standardizedX;
    if (halfNoncentrality == 0)
        return -std::expm1(-halfX);
    return math::MarcumP(1.0, halfNoncentrality, halfX);
}

template <typename RealType>
double RiceRand<RealType>::S(const RealType &x) const
{
    if (x <= 0)
        return 1;
    if (std::isinf(x))
        return 0;
    const double standardizedX = static_cast<double>(x) / sigma;
    const double halfX = 0.5 * standardizedX * standardizedX;
    if (halfNoncentrality == 0)
        return std::exp(-halfX);
    return math::MarcumQ(1.0, halfNoncentrality, halfX);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType RiceRand<RealType>::operator()(G &g) const
{
    const RealType y = static_cast<RealType>(nu) + static_cast<RealType>(sigma) * detail::standardNormal<RealType>(g);
    const RealType z = static_cast<RealType>(sigma) * detail::standardNormal<RealType>(g);
    return std::hypot(y, z);
}

template <typename RealType>
void RiceRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    for (RealType &value : outputData)
        value = (*this)(g);
}

template <typename RealType>
long double RiceRand<RealType>::Mean() const
{
    const long double k = halfNoncentrality;
    const long double h = 0.5L * k;
    /// For large h the scaled Bessel functions e^{-h} I(h) come from their Hankel series: log I(h) - h would
    /// lose h ulps of the logarithm where long double is binary64
    static constexpr long double HANKEL_ARGUMENT = 50;
    long double scaledI0 = 0, ratio = 0;
    if (h > HANKEL_ARGUMENT) {
        const long double logSeries0 = math::detail::logHankelSeries<long double>(0, h, -1);
        scaledI0 = std::exp(logSeries0) / std::sqrt(2 * static_cast<long double>(randlib::consts::pi) * h);
        ratio = std::exp(math::detail::logHankelSeries<long double>(1, h, -1) - logSeries0);
    } else {
        const long double logI0 = math::logBesselI(0.0L, h);
        scaledI0 = std::exp(logI0 - h);
        ratio = (h == 0) ? 0 : std::exp(math::logBesselI(1.0L, h) - logI0);
    }
    return static_cast<long double>(sigma) * std::sqrt(randlib::consts::pi / 2.0L) * scaledI0 * (1 + k + k * ratio);
}

template <typename RealType>
long double RiceRand<RealType>::Variance() const
{
    const long double mean = Mean();
    return 2 * static_cast<long double>(sigma) * sigma + static_cast<long double>(nu) * nu - mean * mean;
}

template <typename RealType>
RealType RiceRand<RealType>::Mode() const
{
    return (nu == 0) ? static_cast<RealType>(sigma) : ContinuousDistribution<RealType>::Mode();
}

template <typename RealType>
std::complex<double> RiceRand<RealType>::CFImpl(double t) const
{
    return this->integrateCharacteristicFunction(t);
}

} // namespace randlib

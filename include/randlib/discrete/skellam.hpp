#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The SkellamRand class <BR>
 * Skellam distribution
 *
 * Notation: X ~ Skellam(μ1, μ2)
 *
 * Related distributions: <BR>
 * If Y ~ Po(μ1) and Z ~ Po(μ2) then Y - Z ~ Skellam(μ1, μ2)
 */
template <typename IntType>
class SkellamRand final : public DiscreteDistribution<IntType>
{
    double mu1 = 1; ///< first rate μ1
    double mu2 = 1; ///< second rate μ2
    double logMu1 = 0; ///< log(μ1)
    double logMu2 = 0; ///< log(μ2)
    double sqrtMu1 = 1; ///< √μ1
    double sqrtMu2 = 1; ///< √μ2

    PoissonRand<IntType> X{}, Y{};

public:
    SkellamRand(double rate1 = 1.0, double rate2 = 1.0);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] IntType MinValue() const override { return std::numeric_limits<IntType>::lowest(); }
    [[nodiscard]] IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetRates(double rate1, double rate2);
    [[nodiscard]] inline double GetFirstRate() const noexcept { return mu1; }
    [[nodiscard]] inline double GetSecondRate() const noexcept { return mu2; }

    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;
    [[nodiscard]] double S(const IntType &k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

protected:
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }

public:
protected:
    void sample(detail::EngineRef g, std::span<IntType> outputData) const override;

public:
    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Median() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/discrete/poisson.hpp>
#include <randlib/continuous/noncentral_chi_squared.hpp>

namespace randlib
{

template <typename IntType>
SkellamRand<IntType>::SkellamRand(double rate1, double rate2)
{
    SetRates(rate1, rate2);
}

template <typename IntType>
std::string SkellamRand<IntType>::Name() const
{
    return "Skellam(" + this->toStringWithPrecision(GetFirstRate()) + ", " +
           this->toStringWithPrecision(GetSecondRate()) + ")";
}

template <typename IntType>
void SkellamRand<IntType>::SetRates(double rate1, double rate2)
{
    if (!(rate1 > 0.0 && rate2 > 0.0))
        throw std::invalid_argument("Skellam distribution: rates should be positive");

    X.SetRate(rate1);
    mu1 = X.GetRate();
    logMu1 = std::log(mu1);
    sqrtMu1 = std::sqrt(mu1);

    Y.SetRate(rate2);
    mu2 = Y.GetRate();
    logMu2 = std::log(mu2);
    sqrtMu2 = std::sqrt(mu2);
}

template <typename IntType>
double SkellamRand<IntType>::P(const IntType &k) const
{
    return std::exp(logP(k));
}

template <typename IntType>
double SkellamRand<IntType>::logP(const IntType &k) const
{
    double y = math::logBesselI(k, 2 * sqrtMu1 * sqrtMu2);
    y += 0.5 * k * (logMu1 - logMu2);
    y -= mu1 + mu2;
    return y;
}

template <typename IntType>
double SkellamRand<IntType>::F(const IntType &k) const
{
    return (k < 0) ? math::MarcumP(-k, mu1, mu2, sqrtMu1, sqrtMu2, logMu1, logMu2)
                   : math::MarcumQ(k + 1, mu2, mu1, sqrtMu2, sqrtMu1, logMu2, logMu1);
}

template <typename IntType>
double SkellamRand<IntType>::S(const IntType &k) const
{
    return (k < 0) ? math::MarcumQ(-k, mu1, mu2, sqrtMu1, sqrtMu2, logMu1, logMu2)
                   : math::MarcumP(k + 1, mu2, mu1, sqrtMu2, sqrtMu1, logMu2, logMu1);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType SkellamRand<IntType>::operator()(G &g) const
{
    return X(g) - Y(g);
}

template <typename IntType>
void SkellamRand<IntType>::sample(detail::EngineRef g, std::span<IntType> outputData) const
{
    X.Sample(g, outputData);
    for (IntType &var : outputData)
        var -= Y(g);
}

template <typename IntType>
long double SkellamRand<IntType>::Mean() const
{
    return mu1 - mu2;
}

template <typename IntType>
long double SkellamRand<IntType>::Variance() const
{
    return mu1 + mu2;
}

template <typename IntType>
IntType SkellamRand<IntType>::Median() const
{
    return DiscreteDistribution<IntType>::quantileImpl(0.5, mu1 - mu2);
}

template <typename IntType>
std::complex<double> SkellamRand<IntType>::CFImpl(double t) const
{
    double cosT = std::cos(t), sinT = std::sin(t);
    double x = (cosT - 1) * (mu1 + mu2);
    double y = sinT * (mu1 - mu2);
    std::complex<double> z(x, y);
    return std::exp(z);
}

template <typename IntType>
IntType SkellamRand<IntType>::Mode() const
{
    return Mean();
}

template <typename IntType>
long double SkellamRand<IntType>::Skewness() const
{
    return (mu1 - mu2) / std::pow(mu1 + mu2, 1.5);
}

template <typename IntType>
long double SkellamRand<IntType>::ExcessKurtosis() const
{
    return 1.0 / (mu1 + mu2);
}

} // namespace randlib

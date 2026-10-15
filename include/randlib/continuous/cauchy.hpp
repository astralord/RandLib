#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/stable.hpp>

namespace randlib
{

/**
 * @brief The CauchyRand class <BR>
 * Cauchy distribution
 *
 * f(x | μ, σ) = σ / [π (σ^2 + (x - μ)^2)]
 *
 * Notation: X ~ Cauchy(μ, σ)
 *
 * Related distributions: <BR>
 * If X ~ Cauchy(0, 1), then μ + σ * X ~ Cauchy(μ, σ) <BR>
 * X ~ S(1, 0, σ, μ) <BR>
 * If X, Y ~ Normal(0, 1), then X / Y ~ Cauchy(0, 1)
 */
template <typename RealType>
class CauchyRand final : public StableDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

public:
    CauchyRand(double location = 0, double scale = 1);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] RealType MinValue() const override { return -INFINITY; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Entropy() const;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template <typename RealType>
CauchyRand<RealType>::CauchyRand(double location, double scale) : StableDistribution<RealType>(1, 0, scale, location)
{
}

template <typename RealType>
std::string CauchyRand<RealType>::Name() const
{
    return "Cauchy(" + this->toStringWithPrecision(this->GetLocation()) + ", " +
           this->toStringWithPrecision(this->GetScale()) + ")";
}

template <typename RealType>
double CauchyRand<RealType>::f(const RealType &x) const
{
    return this->pdfCauchy(x);
}

template <typename RealType>
double CauchyRand<RealType>::F(const RealType &x) const
{
    return this->cdfCauchy(x);
}

template <typename RealType>
double CauchyRand<RealType>::S(const RealType &x) const
{
    return this->cdfCauchyCompl(x);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType CauchyRand<RealType>::operator()(G &g) const
{
    return this->mu + this->gamma * StandardVariate(g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType CauchyRand<RealType>::StandardVariate(G &g)
{
    double x, y;
    do {
        x = 2 * detail::uniform01<RealType>(g) - 1;
        y = 2 * detail::uniform01<RealType>(g) - 1;
    } while (y == 0.0 || x * x + y * y > 1.0);
    return x / y;
}

template <typename RealType>
std::complex<double> CauchyRand<RealType>::CFImpl(double t) const
{
    return this->cfCauchy(t);
}

template <typename RealType>
RealType CauchyRand<RealType>::quantileImpl(double p) const
{
    return this->quantileCauchy(p);
}

template <typename RealType>
RealType CauchyRand<RealType>::quantileImpl1m(double p) const
{
    return this->quantileCauchy1m(p);
}

template <typename RealType>
long double CauchyRand<RealType>::Entropy() const
{
    return 2 * randlib::consts::ln2 + this->logGamma + randlib::consts::lnpi;
}

} // namespace randlib

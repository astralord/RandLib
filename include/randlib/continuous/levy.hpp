#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/stable.hpp>

namespace randlib
{

/**
 * @brief The LevyRand class <BR>
 * Levy distribution
 *
 * f(x | μ, σ) = ((σ exp(σ / (μ - x)) / (2 π (x - μ)^3))^(1/2)
 *
 * Notation: X ~ Levy(μ, σ)
 *
 * Related distributions: <BR>
 * If X ~ Levy(0, 1), then μ + σ * X ~ Levy(μ, σ) <BR>
 * X ~ S(0.5, 1, σ, μ) <BR>
 * If Y ~ Normal(0, 1), then 1 / X^2 ~ Levy(0, 1) <BR>
 * If X ~ Levy(0, σ), then X ~ Inv-Γ(1/2, σ/2)
 */
template <typename RealType>
class LevyRand final : public StableDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

public:
    LevyRand(double location = 0, double scale = 1);
    [[nodiscard]] std::string Name() const override;

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
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
    /**
     * Fit scale using maximum-likelihoood estimator
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitScale(std::span<const RealType> sample);
};

} // namespace randlib

#include <randlib/continuous/normal.hpp>

namespace randlib
{

template <typename RealType>
LevyRand<RealType>::LevyRand(double location, double scale) : StableDistribution<RealType>(0.5, 1, scale, location)
{
}

template <typename RealType>
std::string LevyRand<RealType>::Name() const
{
    return "Levy(" + this->toStringWithPrecision(this->GetLocation()) + ", " +
           this->toStringWithPrecision(this->GetScale()) + ")";
}

template <typename RealType>
double LevyRand<RealType>::f(const RealType &x) const
{
    return this->pdfLevy(x);
}

template <typename RealType>
double LevyRand<RealType>::logf(const RealType &x) const
{
    return this->logpdfLevy(x);
}

template <typename RealType>
double LevyRand<RealType>::F(const RealType &x) const
{
    return this->cdfLevy(x);
}

template <typename RealType>
double LevyRand<RealType>::S(const RealType &x) const
{
    return this->cdfLevyCompl(x);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType LevyRand<RealType>::operator()(G &g) const
{
    RealType rv = detail::standardNormal<RealType>(g);
    rv *= rv;
    rv = this->gamma / rv;
    return this->mu + rv;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType LevyRand<RealType>::StandardVariate(G &g)
{
    RealType rv = detail::standardNormal<RealType>(g);
    return 1.0 / (rv * rv);
}

template <typename RealType>
RealType LevyRand<RealType>::quantileImpl(double p) const
{
    return this->quantileLevy(p);
}

template <typename RealType>
RealType LevyRand<RealType>::quantileImpl1m(double p) const
{
    return this->quantileLevy1m(p);
}

template <typename RealType>
std::complex<double> LevyRand<RealType>::CFImpl(double t) const
{
    return this->cfLevy(t);
}

template <typename RealType>
void LevyRand<RealType>::FitScale(std::span<const RealType> sample)
{
    /// Sanity check
    if (!this->allElementsAreNotSmallerThan(this->mu, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->mu));
    long double invSum = 0.0;
    for (RealType var : sample)
        invSum += 1.0 / (var - this->mu);
    invSum = 1.0 / invSum;
    this->SetScale(sample.size() * invSum);
}

} // namespace randlib

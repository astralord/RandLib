#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/circular.hpp>

namespace randlib
{

/**
 * @brief The WrappedExponentialRand class <BR>
 * Wrapped Exponential distribution
 *
 * Notation: X ~ WE(λ)
 */
template <typename RealType>
class WrappedExponentialRand final : public CircularDistribution<RealType>
{
    double lambda = 1; ///< rate λ
    double logLambda = 0; ///< ln(λ)
    double scaledLambda = 2 * randlib::consts::pi; ///< 2πλ
    double pdfCoef = 0, logpdfCoef = 0, expmScaledLambda = 0;

public:
    WrappedExponentialRand(double rate);

    [[nodiscard]] std::string Name() const override;

    void SetRate(double rate);
    [[nodiscard]] inline double GetRate() const noexcept { return lambda; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double CircularMean() const override;
    [[nodiscard]] long double CircularVariance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;

protected:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template < typename RealType >
WrappedExponentialRand<RealType>::WrappedExponentialRand(double rate) : CircularDistribution<RealType>(randlib::consts::pi)
{
    SetRate(rate);
}

template < typename RealType >
std::string WrappedExponentialRand<RealType>::Name() const
{
  return "Wrapped Exponential(" + this->toStringWithPrecision(GetRate()) + ")";
}

template < typename RealType >
void WrappedExponentialRand<RealType>::SetRate(double rate)
{
    if (!(rate > 0.0))
        throw std::invalid_argument("Wrapped Exponential distribution: rate parameter should be positive");
    lambda = rate;
    logLambda = std::log(lambda);
    scaledLambda = 2 * randlib::consts::pi * lambda;
    expmScaledLambda = std::exp(-scaledLambda);
    pdfCoef = -std::expm1l(-scaledLambda);
    logpdfCoef = math::log1mexp(-scaledLambda);
}

template < typename RealType >
double WrappedExponentialRand<RealType>::f(const RealType &x) const
{
  return (x < this->MinValue() || x > this->MaxValue()) ? 0.0 : std::exp(logf(x));
}

template < typename RealType >
double WrappedExponentialRand<RealType>::logf(const RealType &x) const
{
  return (x < this->MinValue() || x > this->MaxValue()) ? - INFINITY : logLambda - lambda * x - logpdfCoef;
}

template < typename RealType >
double WrappedExponentialRand<RealType>::F(const RealType &x) const
{
    if (x <= 0.0)
      return 0.0;
  return (x < this->MaxValue()) ? std::exp(math::log1mexp(-lambda * x) - logpdfCoef) : 1.0;
}

template < typename RealType >
double WrappedExponentialRand<RealType>::S(const RealType &x) const
{
    if (x <= 0.0)
      return 1.0;
    if (x >= this->MaxValue())
      return 0.0;
    double y = std::expm1l(scaledLambda - lambda * x);
    y /= pdfCoef;
  return expmScaledLambda * y;
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType WrappedExponentialRand<RealType>::operator()(G &g) const
{
  return quantileImpl(detail::uniform01<RealType>(g));
}

template < typename RealType >
long double WrappedExponentialRand<RealType>::CircularMean() const
{
  return randlib::consts::half_pi - math::atan(lambda);
}

template < typename RealType >
long double WrappedExponentialRand<RealType>::CircularVariance() const
{
  return 1.0 - lambda / std::sqrt(1.0 + lambda * lambda);
}

template < typename RealType >
RealType WrappedExponentialRand<RealType>::Median() const
{
  return (randlib::consts::ln2 - math::softplus(-scaledLambda)) / lambda;
}

template < typename RealType >
RealType WrappedExponentialRand<RealType>::Mode() const
{
  return 0.0;
}

template < typename RealType >
RealType WrappedExponentialRand<RealType>::quantileImpl(double p) const
{
  return -std::log1pl(-p * pdfCoef) / lambda;
}

template < typename RealType >
RealType WrappedExponentialRand<RealType>::quantileImpl1m(double p) const
{
  return -std::log(expmScaledLambda + p * pdfCoef) / lambda;
}

template < typename RealType >
std::complex<double> WrappedExponentialRand<RealType>::CFImpl(double t) const
{
    /// φ(t) = λ (1 - e^{2π(it - λ)}) / ((λ - it)(1 - e^{-2πλ})), which reduces to λ / (λ - it) for integer t
    std::complex<double> lambdaMinusIt(lambda, -t);
    std::complex<double> wrap = 1.0 - std::exp(-lambdaMinusIt * static_cast<double>(2 * randlib::consts::pi));
  return lambda * wrap / (lambdaMinusIt * pdfCoef);
}

} // namespace randlib

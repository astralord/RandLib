#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The SechRand class <BR>
 * Hyperbolic secant distribution
 *
 * Notation: X ~ Sech
 */
template <typename RealType>
class SechRand final : public ContinuousDistribution<RealType>
{
public:
    SechRand();

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] RealType MinValue() const override { return -INFINITY; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;
public:
    [[nodiscard]] long double Entropy() const;
};

} // namespace randlib

#include <randlib/continuous/cauchy.hpp>

namespace randlib
{

template < typename RealType >
SechRand<RealType>::SechRand()
{
}

template < typename RealType >
std::string SechRand<RealType>::Name() const
{
  return "Hyperbolic secant";
}

template < typename RealType >
double SechRand<RealType>::f(const RealType & x) const
{
  return 0.5 / std::cosh(randlib::consts::half_pi * x);
}

template < typename RealType >
double SechRand<RealType>::logf(const RealType & x) const
{
  return randlib::consts::half_pi * x - math::softplus(randlib::consts::pi * x);
}

template < typename RealType >
double SechRand<RealType>::F(const RealType & x) const
{
    double y = std::exp(randlib::consts::half_pi * x);
  return randlib::consts::two_inv_pi * math::atan(y);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType SechRand<RealType>::operator()(G &g) const
{
    RealType y = std::fabs(CauchyRand<RealType>::StandardVariate(g));
  return randlib::consts::two_inv_pi * std::log(y);
}

template < typename RealType >
long double SechRand<RealType>::Mean() const
{
  return 0.0;
}

template < typename RealType >
long double SechRand<RealType>::Variance() const
{
  return 1.0;
}

template < typename RealType >
std::complex<double> SechRand<RealType>::CFImpl(double t) const
{
  return 1.0 / std::cosh(t);
}

template < typename RealType >
RealType SechRand<RealType>::quantileImpl(double p) const
{
    RealType x = randlib::consts::half_pi * p;
    x = std::tan(x);
    x = std::log(x);
  return randlib::consts::two_inv_pi * x;
}

template < typename RealType >
RealType SechRand<RealType>::quantileImpl1m(double p) const
{
    RealType x = randlib::consts::half_pi * p;
    x = std::tan(x);
    x = -std::log(x);
  return randlib::consts::two_inv_pi * x;
}

template < typename RealType >
RealType SechRand<RealType>::Median() const
{
  return 0.0;
}

template < typename RealType >
RealType SechRand<RealType>::Mode() const
{
  return 0.0;
}

template < typename RealType >
long double SechRand<RealType>::Skewness() const
{
  return 0.0;
}

template < typename RealType >
long double SechRand<RealType>::ExcessKurtosis() const
{
  return 2.0;
}

template < typename RealType >
long double SechRand<RealType>::Entropy() const
{
  return 2.0 * randlib::consts::two_inv_pi * randlib::consts::catalan;
}

} // namespace randlib

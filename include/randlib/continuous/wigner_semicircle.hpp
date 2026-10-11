#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/beta.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The WignerSemicircleRand class <BR>
 * Wigner-Semicircle distribution
 *
 * Notation: X ~ Wigner-Sc(R)
 *
 * Related distributions:
 * If Y ~ Beta(1.5, 1.5), then R * (2Y - 1) ~ Wigner-Sc(R)
 */
template <typename RealType>
class WignerSemicircleRand final : public ContinuousDistribution<RealType>
{
    RealType R = 1; ///< radius
    double RSq = 1; ///< R^2
    double logRSq = 0; /// log(R^2)
    BetaRand<RealType> X{1.5, 1.5};
    
    double lowerTail(double x) const;

public:
    explicit WignerSemicircleRand(double radius);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return -R; }
    [[nodiscard]] RealType MaxValue() const override { return R; }

    void SetRadius(double radius);
    [[nodiscard]] inline double GetRadius() const noexcept { return R; }
    
public:
    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;
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
    
    [[nodiscard]] long double Entropy() const;
};

} // namespace randlib

namespace randlib
{

template < typename RealType >
WignerSemicircleRand<RealType>::WignerSemicircleRand(double radius)
{
    SetRadius(radius);
}

template < typename RealType >
std::string WignerSemicircleRand<RealType>::Name() const
{
  return "Wigner Semicircle(" + this->toStringWithPrecision(GetRadius()) + ")";
}

template < typename RealType >
void WignerSemicircleRand<RealType>::SetRadius(double radius)
{
    if (!(radius > 0.0))
        throw std::invalid_argument("Wigner-Semicircle distribution: radius should be positive");
    R = radius;
    RSq = R * R;
    logRSq = std::log(RSq);
}

template < typename RealType >
double WignerSemicircleRand<RealType>::f(const RealType & x) const
{
    double xSq = x * x;
    if (xSq >= RSq)
      return 0.0;
    double y = RSq - xSq;
    y = std::sqrt(y);
    y *= randlib::consts::inv_pi / RSq;
  return 2 * y;
}

template < typename RealType >
double WignerSemicircleRand<RealType>::logf(const RealType & x) const
{
    double xSq = x * x;
    if (xSq >= RSq)
      return -INFINITY;
  return randlib::consts::ln2 + 0.5 * std::log(RSq - xSq) - randlib::consts::lnpi - logRSq;
}

template < typename RealType >
double WignerSemicircleRand<RealType>::lowerTail(double x) const
{
    /// F(x) for x ≤ 0. With x = R sin θ and φ = θ + π/2, F = (2φ - sin 2φ) / (2π); taking
    /// φ = 2 asin(√(δ/2)) with δ = 1 + x/R avoids the cancellation of the textbook formula near -R
    double delta = (R + x) / R;
    double phi = 2 * std::asin(std::sqrt(0.5 * delta));
    double y = 2 * phi, ySq = y * y;
    double yMinusSinY = (y < 0.5) ? y * ySq / 6 * (1 - ySq / 20 * (1 - ySq / 42 * (1 - ySq / 72 * (1 - ySq / 110))))
                                  : y - std::sin(y);
  return yMinusSinY / (2 * randlib::consts::pi);
}

template < typename RealType >
double WignerSemicircleRand<RealType>::F(const RealType &x) const
{
    if (x <= -R)
      return 0.0;
    if (x >= R)
      return 1.0;
  return (x <= 0) ? lowerTail(x) : 1.0 - lowerTail(-x);
}

template < typename RealType >
double WignerSemicircleRand<RealType>::S(const RealType &x) const
{
    if (x <= -R)
      return 1.0;
    if (x >= R)
      return 0.0;
  return (x >= 0) ? lowerTail(-x) : 1.0 - lowerTail(x);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType WignerSemicircleRand<RealType>::operator()(G &g) const
{
    RealType x = X(g);
    x += x - 1;
  return R * x;
}



template < typename RealType >
long double WignerSemicircleRand<RealType>::Mean() const
{
  return 0.0;
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::Variance() const
{
  return 0.25 * RSq;
}

template < typename RealType >
RealType WignerSemicircleRand<RealType>::Median() const
{
  return 0.0;
}

template < typename RealType >
RealType WignerSemicircleRand<RealType>::Mode() const
{
  return 0.0;
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::Skewness() const
{
  return 0.0;
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::ExcessKurtosis() const
{
  return -1.0;
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::Entropy() const
{
  return randlib::consts::lnpi + 0.5 * logRSq - 0.5;
}

} // namespace randlib

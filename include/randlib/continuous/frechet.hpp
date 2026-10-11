#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The FrechetRand class <BR>
 * Frechet distribution
 *
 * Notation: X ~ Frechet(α, s, m)
 */
template <typename RealType>
class FrechetRand final : public ContinuousDistribution<RealType>
{
    double alpha = 1; ///< shape α
    double s = 1; ///< scale
    double m = 0; ///< location
    double alphaInv = 1; ///< 1/α
    double pdfCoef = 0; ///< log(α/s)

public:
    FrechetRand(double shape = 1, double scale = 1, double location = 0);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return m; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double scale, double location);
    [[nodiscard]] inline double GetShape() const noexcept { return alpha; }
    [[nodiscard]] inline double GetScale() const noexcept { return s; }
    [[nodiscard]] inline double GetLocation() const noexcept { return m; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;
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

public:
    [[nodiscard]] long double Entropy() const;
};

} // namespace randlib

#include <randlib/continuous/exponential.hpp>

namespace randlib
{

template < typename RealType >
FrechetRand<RealType>::FrechetRand(double shape, double scale, double location)
{
    SetParameters(shape, scale, location);
}

template < typename RealType >
std::string FrechetRand<RealType>::Name() const
{
  return "Frechet(" + this->toStringWithPrecision(GetShape()) + ", "
                      + this->toStringWithPrecision(GetScale()) + ", "
                      + this->toStringWithPrecision(GetLocation()) + ")";
}

template < typename RealType >
void FrechetRand<RealType>::SetParameters(double shape, double scale, double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Frechet distribution: location should be finite");
    if (!(shape > 0.0))
        throw std::invalid_argument("Frechet distribution: shape should be positive");
    if (!(scale > 0.0))
        throw std::invalid_argument("Frechet distribution: scale should be positive");
    alpha = shape;
    alphaInv = 1.0 / alpha;
    s = scale;
    pdfCoef = std::log(alpha / s);
    m = location;
}

template < typename RealType >
double FrechetRand<RealType>::f(const RealType &x) const
{
  return (x <= m) ? 0.0 : std::exp(logf(x));
}

template < typename RealType >
double FrechetRand<RealType>::logf(const RealType &x) const
{
    if (x <= m)
      return -INFINITY;
    double xAdj = (x - m) / s;
    double logxAdj = std::log(xAdj);
    double a = alpha * logxAdj;
    double expA = std::exp(-a);
  return pdfCoef - a - expA - logxAdj;
}

template < typename RealType >
double FrechetRand<RealType>::F(const RealType & x) const
{
    if (x <= m)
      return 0.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
  return std::exp(-xPow);
}

template < typename RealType >
double FrechetRand<RealType>::S(const RealType & x) const
{
    if (x <= m)
      return 1.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
  return -std::expm1l(-xPow);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType FrechetRand<RealType>::operator()(G &g) const
{
  return m + s / std::pow(detail::standardExponential<RealType>(g), alphaInv);
}

template < typename RealType >
long double FrechetRand<RealType>::Mean() const
{
    if (alpha <= 1.0)
      return INFINITY;
  return m + s * std::tgammal(1.0 - alphaInv);
}

template < typename RealType >
long double FrechetRand<RealType>::Variance() const
{
    if (alpha <= 2.0)
      return INFINITY;
    long double var = std::tgammal(1.0 - alphaInv);
    var *= var;
    var = std::tgammal(1.0 - 2 * alphaInv) - var;
  return s * s * var;
}

template < typename RealType >
RealType FrechetRand<RealType>::quantileImpl(double p) const
{
    RealType y = -std::log(p);
    y = s / std::pow(y, alphaInv);
  return y + m;
}

template < typename RealType >
RealType FrechetRand<RealType>::quantileImpl1m(double p) const
{
    RealType y = -std::log1pl(-p);
    y = s / std::pow(y, alphaInv);
  return y + m;
}

template < typename RealType >
RealType FrechetRand<RealType>::Median() const
{
  return m + s / std::pow(randlib::consts::ln2, alphaInv);
}

template < typename RealType >
RealType FrechetRand<RealType>::Mode() const
{
    RealType y = alpha / (1.0 + alpha);
    y = std::pow(y, alphaInv);
  return m + s * y;
}

template < typename RealType >
long double FrechetRand<RealType>::Skewness() const
{
    if (alpha <= 3.0)
      return INFINITY;
    long double x = std::tgammal(1.0 - alphaInv);
    long double y = std::tgammal(1.0 - 2.0 * alphaInv);
    long double z = std::tgammal(1.0 - 3.0 * alphaInv);
    long double numerator = 2 * x * x - 3 * y;
    numerator *= x;
    numerator += z;
    long double denominator = y - x * x;
    denominator = std::pow(denominator, 1.5);
  return numerator / denominator;
}

template < typename RealType >
long double FrechetRand<RealType>::ExcessKurtosis() const
{
    if (alpha <= 4.0)
      return INFINITY;
    long double x = std::tgammal(1.0 - alphaInv);
    long double y = std::tgammal(1.0 - 2.0 * alphaInv);
    long double z = std::tgammal(1.0 - 3.0 * alphaInv);
    long double w = std::tgammal(1.0 - 4.0 * alphaInv);
    long double numerator = w - 4 * z * x + 3 * y * y;
    long double denominator = y - x * x;
    denominator *= denominator;
  return numerator / denominator - 6.0;
}

template < typename RealType >
long double FrechetRand<RealType>::Entropy() const
{
  return 1.0 + randlib::consts::egamma * (1.0 + alphaInv) + std::log(s / alpha);
}

} // namespace randlib

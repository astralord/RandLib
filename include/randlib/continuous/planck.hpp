#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/gamma.hpp>
#include <randlib/discrete/zeta.hpp>

namespace randlib
{

/**
 * @brief The PlanckRand class <BR>
 * Planck distribution
 *
 * f(x | a, b) = g(a + 1) * (x ^ a) / (exp(b * x) - 1)
 * where g(y) = b ^ y / (Γ(y) * ζ(y))
 *
 * Notation: X ~ Planck(a, b)
 *
 * Related distributions: <BR>
 * If G ~ Gamma(a + 1, b) and Z ~ Zeta(a + 1), then G / Z ~ Planck(a, b)
 */
template <typename RealType>
class PlanckRand final : public ContinuousDistribution<RealType>
{
    double a = 1; ///< shape
    double b = 1; ///< scale
    double pdfCoef = randlib::consts::ln2 + randlib::consts::ln3 - 2 * randlib::consts::lnpi; ///< coefficient for faster pdf calculations

    ZetaRand<long long int> Z{};
    GammaRand<RealType> gammaVariable{2};

public:
    PlanckRand(double shape = 1, double scale = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double scale);
    [[nodiscard]] inline double GetShape() const noexcept { return a; }
    [[nodiscard]] inline double GetScale() const noexcept { return b; }

private:
    /**
     * @fn h
     * @param t
     * @return leveled pdf for F and CF calculations
     */
    double h(double t) const;

public:
    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:
    protected:
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double SecondMoment() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double ThirdMoment() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double FourthMoment() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

namespace randlib
{

template < typename RealType >
PlanckRand<RealType>::PlanckRand(double shape, double scale)
{
    SetParameters(shape, scale);
}

template < typename RealType >
std::string PlanckRand<RealType>::Name() const
{
  return "Planck(" + this->toStringWithPrecision(GetShape()) + ", " + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void PlanckRand<RealType>::SetParameters(double shape, double scale)
{
    if (!(shape > 0.0))
        throw std::invalid_argument("Planck distribution: shape should be positive, but it's equal to " + std::to_string(shape));
    if (!(scale > 0.0))
        throw std::invalid_argument("Planck distribution: scale should be positive, but it's equal to " + std::to_string(scale));
    a = shape;
    b = scale;

    double ap1 = a + 1;
    Z.SetExponent(ap1);
    gammaVariable.SetParameters(ap1, b);

    pdfCoef = -Z.GetLogZetaFunction();
    pdfCoef += ap1 * std::log(b);
    pdfCoef -= gammaVariable.GetLogGammaShape();
}

template < typename RealType >
double PlanckRand<RealType>::h(double t) const
{
    if (t <= 0)
      return 0.0;
    double y = pdfCoef + a * std::log(t);
    double expY = std::exp(y);
    double bt = b * t;
    double z = 1.0 / std::expm1l(bt) - 1.0 / bt;
  return expY * z;
}

template < typename RealType >
double PlanckRand<RealType>::f(const RealType &x) const
{
    if (x < 0)
      return 0;
    if (x == 0)
    {
        if (a > 1)
          return 0.0;
        if (a == 1)
          return std::exp(pdfCoef) / b;
      return INFINITY;
    }
    double y = pdfCoef + a * std::log(x);
  return std::exp(y) / std::expm1l(b * x);
}

template < typename RealType >
double PlanckRand<RealType>::logf(const RealType &x) const
{
    if (x < 0.0)
      return -INFINITY;
    if (x == 0)
    {
        if (a > 1)
          return -INFINITY;
        if (a == 1)
          return pdfCoef - gammaVariable.GetLogRate();
      return INFINITY;
    }
    double y = pdfCoef + a * std::log(x);
  return y - math::logexpm1l(b * x);
}

template < typename RealType >
double PlanckRand<RealType>::F(const RealType &x) const
{
    if (x <= 0)
      return 0.0;

    if (a >= 1) {
      return math::integral([this] (double t)
        {
            return f(t);
        },
        0, x);
    }

    /// split F(x) by two integrals
    double aux = pdfCoef + a * std::log(x);
    double integral1 = std::exp(aux) / (b * a);
    double integral2 = math::integral([this] (double t)
    {
      return h(t);
    },
    0, x);
  return integral1 + integral2;
}



template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType PlanckRand<RealType>::operator()(G &g) const
{
  return gammaVariable(g) / Z(g);
}

template < typename RealType >
void PlanckRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    gammaVariable.Sample(g, outputData);
    for (RealType & var : outputData)
        var /= Z(g);
}

template < typename RealType >
long double PlanckRand<RealType>::Mean() const
{
    double y = (a + 1) / b;
    y *= math::riemannZeta(a + 2);
  return y / Z.GetZetaFunction();
}

template < typename RealType >
long double PlanckRand<RealType>::SecondMoment() const
{
    long double secondMoment = (a + 1) * (a + 2);
    secondMoment /= (b * b);
    secondMoment *= math::riemannZeta(a + 3);
    secondMoment /= Z.GetZetaFunction();
  return secondMoment;
}

template < typename RealType >
long double PlanckRand<RealType>::Variance() const
{
    long double mean = Mean();
    return SecondMoment() - mean * mean;
}

template < typename RealType >
RealType PlanckRand<RealType>::Mode() const
{
    if (a <= 1)
      return 0.0;
    double y = -a * std::exp(-a);
    y = math::W0Lambert(y);
  return (y + a) / b;
}

template < typename RealType >
long double PlanckRand<RealType>::ThirdMoment() const
{
    long double thirdMoment = (a + 3) * (a + 2) * (a + 1);
    thirdMoment /= (b * b * b);
    thirdMoment *= math::riemannZeta(a + 4);
    thirdMoment /= Z.GetZetaFunction();
  return thirdMoment;
}

template < typename RealType >
long double PlanckRand<RealType>::Skewness() const
{
    long double mean = Mean();
    long double secondMoment = SecondMoment();
    long double thirdMoment = ThirdMoment();
    long double meanSq = mean * mean;
    long double variance = secondMoment - meanSq;
    long double numerator = thirdMoment - 3 * mean * variance - mean * meanSq;
    long double denominator = std::pow(variance, 1.5);
  return numerator / denominator;
}

template < typename RealType >
long double PlanckRand<RealType>::FourthMoment() const
{
    long double fourthMoment = (a + 4) * (a + 3) * (a + 2) * (a + 1);
    long double bSq = b * b;
    fourthMoment /= (bSq * bSq);
    fourthMoment *= math::riemannZeta(a + 5);
    fourthMoment /= Z.GetZetaFunction();
  return fourthMoment;
}

template < typename RealType >
long double PlanckRand<RealType>::ExcessKurtosis() const
{
    long double mean = Mean();
    long double secondMoment = SecondMoment();
    long double thirdMoment = ThirdMoment();
    long double fourthMoment = FourthMoment();
    long double meanSq = mean * mean;
    long double variance = secondMoment - meanSq;
    long double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    long double denominator = variance * variance;
  return numerator / denominator - 3.0;
}

template < typename RealType >
std::complex<double> PlanckRand<RealType>::CFImpl(double t) const
{
    if (a >= 1)
      return ContinuousDistribution<RealType>::CFImpl(t);

    /// We have singularity point at 0 for real part,
    /// so we split the integral in two intervals:
    /// First one from 0 to 1, for which we integrate
    /// numerically leveled pdf and add known solution for level.
    /// Second one from 1 to infinity, for which we use
    /// simple expected value for the rest of the function
    double re1 = math::integral([this, t] (double x)
    {
      return std::cos(t * x) * h(x);
    },
    0.0, 1.0);

    double re2 = this->ExpectedValue([t] (double x)
    {
      return std::cos(t * x);
    },
    1.0, INFINITY);

    double re3 = t * math::integral([this, t] (double x)
    {
        if (x <= 0.0)
          return 0.0;
      return std::sin(t * x) * std::pow(x, a);
    },
    0.0, 1.0);

    re3 += std::cos(t);
    re3 *= std::exp(pdfCoef) / (b * a);

    double im = this->ExpectedValue([t] (double x)
    {
      return std::sin(t * x);
    },
    0.0, INFINITY);

    double re = re1 + re2 + re3;
  return std::complex<double>(re, im);
}

} // namespace randlib

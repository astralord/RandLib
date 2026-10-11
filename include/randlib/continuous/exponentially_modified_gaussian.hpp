#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/continuous/normal.hpp>

namespace randlib
{

/**
 * @brief The ExponentiallyModifiedGaussianRand class <BR>
 * Exponentially-modified Gaussian distribution
 *
 * Notation: X ~ EMG(μ, σ, β)
 *
 * Related distributions: <BR>
 * X = Y + Z, where Y ~ Normal(μ, σ) and Z ~ Exp(β)
 */
template <typename RealType>
class ExponentiallyModifiedGaussianRand final : public ContinuousDistribution<RealType>
{
    NormalRand<RealType> X{};
    ExponentialRand<RealType> Y{};

    double a = 1; ///< μ + βσ^2
    double b = randlib::consts::inv_sqrt2; ///< 1 / (√2 * σ)
    double c = 0.5; ///< μ + βσ^2 / 2
    double v = 1; /// βσ

public:
    explicit ExponentiallyModifiedGaussianRand(double location = 0, double variance = 1, double rate = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] RealType MinValue() const override { return -INFINITY; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double location, double variance, double rate);
    [[nodiscard]] inline double GetLocation() const { return X.GetLocation(); }
    [[nodiscard]] inline double GetScale() const { return X.GetScale(); }
    [[nodiscard]] inline double GetRate() const { return Y.GetRate(); }

private:
    DoublePair faux(const RealType &x) const;

public:
    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

namespace randlib
{

template < typename RealType >
ExponentiallyModifiedGaussianRand<RealType>::ExponentiallyModifiedGaussianRand(double location, double variance, double rate)
{
    SetParameters(location, variance, rate);
}

template < typename RealType >
std::string ExponentiallyModifiedGaussianRand<RealType>::Name() const
{
  return "Exponentially modified Gaussian(" + this->toStringWithPrecision(GetLocation()) + ", "
                                              + this->toStringWithPrecision(X.Variance()) + ", "
                                              + this->toStringWithPrecision(GetRate()) + ")";
}

template < typename RealType >
void ExponentiallyModifiedGaussianRand<RealType>::SetParameters(double location, double variance, double rate)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Exponentially modified Gaussian distribution: location should be finite");
    if (!(variance > 0))
        throw std::invalid_argument("Exponentially modified Gaussian distribution: variance should be positive");
    if (!(rate > 0))
        throw std::invalid_argument("Exponentially modified Gaussian distribution: rate should be positive");

    X.SetLocation(location);
    X.SetVariance(variance);
    Y.SetRate(rate);

    double mu = X.GetLocation();
    double sigma = X.GetScale();
    double beta = Y.GetRate();
    double var = sigma * sigma;
    a = 0.5 * beta * var;
    c = mu + a;
    a += c;
    b = randlib::consts::inv_sqrt2 / sigma;
    v = beta * sigma;
}

template < typename RealType >
DoublePair ExponentiallyModifiedGaussianRand<RealType>::faux(const RealType &x) const
{
    double lambda = Y.GetRate();
    double y = a - x;
    y *= b;
    y = std::erfc(y);
    y *= 0.5 * lambda;
    double exponent = c - x;
    exponent *= lambda;
  return std::make_pair(y, exponent);
}

template < typename RealType >
double ExponentiallyModifiedGaussianRand<RealType>::f(const RealType & x) const
{
    auto [y, exponent] = faux(x);
  return y * std::exp(exponent);
}

template < typename RealType >
double ExponentiallyModifiedGaussianRand<RealType>::logf(const RealType & x) const
{
    auto [y, exponent] = faux(x);
  return std::log(y) + exponent;
}

template < typename RealType >
double ExponentiallyModifiedGaussianRand<RealType>::F(const RealType & x) const
{
    double u = Y.GetRate() * (x - X.GetLocation());
    double y = X.F(x);
    double exponent = -u + 0.5 * v * v;
    exponent = std::exp(exponent);
    exponent *= X.F(x - v * X.GetScale());
  return y - exponent;
}

template < typename RealType >
double ExponentiallyModifiedGaussianRand<RealType>::S(const RealType & x) const
{
    double u = Y.GetRate() * (x - X.GetLocation());
    double y = X.S(x);
    double exponent = -u + 0.5 * v * v;
    exponent = std::exp(exponent);
    exponent *= X.F(x - v * X.GetScale());
  return y + exponent;
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType ExponentiallyModifiedGaussianRand<RealType>::operator()(G &g) const
{
  return X(g) + Y(g);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType ExponentiallyModifiedGaussianRand<RealType>::StandardVariate(G &g)
{
  return detail::standardNormal<RealType>(g) + detail::standardExponential<RealType>(g);
}



template < typename RealType >
long double ExponentiallyModifiedGaussianRand<RealType>::Mean() const
{
  return X.Mean() + Y.Mean();
}

template < typename RealType >
long double ExponentiallyModifiedGaussianRand<RealType>::Variance() const
{
  return X.Variance() + Y.Variance();
}

template < typename RealType >
std::complex<double> ExponentiallyModifiedGaussianRand<RealType>::CFImpl(double t) const
{
  return X.CF(t) * Y.CF(t);
}

template < typename RealType >
long double ExponentiallyModifiedGaussianRand<RealType>::Skewness() const
{
    long double sigma = X.GetScale();
    long double lambda = Y.GetRate();
    long double tmp = 1.0 / (sigma * lambda);
    long double tmpSq = tmp * tmp;
    long double y = 1.0 + tmpSq;
    y = y * y * y;
    y = std::sqrt(y);
    y = tmpSq * tmp / y;
  return y + y;
}

template < typename RealType >
long double ExponentiallyModifiedGaussianRand<RealType>::ExcessKurtosis() const
{
    long double sigma = X.GetScale();
    long double lambda = Y.GetRate();
    long double tmp = 1.0 / (sigma * lambda);
    tmp *= tmp;
    long double numerator = 1.0 + 2.0 * tmp + 3.0 * tmp * tmp;
    long double denominator = 1.0 + tmp;
    denominator *= denominator;
    long double y = numerator / denominator - 1.0;
  return 3.0 * y;
}

} // namespace randlib

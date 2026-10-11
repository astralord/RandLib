#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The GumbelRand class <BR>
 * Gumbel distribution
 *
 * Notation: X ~ Gumbel(μ, β)
 *
 * Related distributions: <BR>
 * exp(-(X - μ) / β) ~ Exp(1)
 */
template <typename RealType>
class GumbelRand final : public ContinuousDistribution<RealType>
{
    double mu = 0; ///< location μ
    double beta = 1; ///< scale β
    double logBeta = 0; ///< log(β)
public:
    GumbelRand(double location = 0, double scale = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] RealType MinValue() const override { return -INFINITY; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    [[nodiscard]] inline double GetLocation() const noexcept { return mu; }
    [[nodiscard]] inline double GetScale() const noexcept { return beta; }

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
GumbelRand<RealType>::GumbelRand(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

template < typename RealType >
std::string GumbelRand<RealType>::Name() const
{
  return "Gumbel(" + this->toStringWithPrecision(GetLocation()) + ", " + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void GumbelRand<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Gumbel distribution: location should be finite");
    mu = location;
}

template < typename RealType >
void GumbelRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Gumbel distribution: scale should be positive");
    beta = scale;
    logBeta = std::log(beta);
}

template < typename RealType >
double GumbelRand<RealType>::f(const RealType & x) const
{
  return std::exp(logf(x));
}

template < typename RealType >
double GumbelRand<RealType>::logf(const RealType & x) const
{
    double z = (mu - x) / beta;
    double y = std::exp(z);
  return z - y - logBeta;
}

template < typename RealType >
double GumbelRand<RealType>::F(const RealType & x) const
{
    double y = (mu - x) / beta;
    y = std::exp(y);
  return std::exp(-y);
}

template < typename RealType >
double GumbelRand<RealType>::S(const RealType &x) const
{
    double y = (mu - x) / beta;
    y = std::exp(y);
  return -std::expm1l(-y);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType GumbelRand<RealType>::operator()(G &g) const
{
  return mu + beta * GumbelRand<RealType>::StandardVariate(g);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType GumbelRand<RealType>::StandardVariate(G &g)
{
    RealType w = detail::standardExponential<RealType>(g);
  return -std::log(w);
}

template < typename RealType >
long double GumbelRand<RealType>::Mean() const
{
  return mu + beta * randlib::consts::egamma;
}

template < typename RealType >
long double GumbelRand<RealType>::Variance() const
{
    double v = randlib::consts::pi * beta;
  return v * v / 6;
}

template < typename RealType >
RealType GumbelRand<RealType>::quantileImpl(double p) const
{
  return mu - beta * std::log(-std::log(p));
}

template < typename RealType >
RealType GumbelRand<RealType>::quantileImpl1m(double p) const
{
  return mu - beta * std::log(-std::log1pl(-p));
}

template < typename RealType >
RealType GumbelRand<RealType>::Median() const
{
    static const double M_LN_LN2 = std::log(randlib::consts::ln2);
  return mu - beta * M_LN_LN2;
}

template < typename RealType >
RealType GumbelRand<RealType>::Mode() const
{
  return mu;
}

template < typename RealType >
long double GumbelRand<RealType>::Skewness() const
{
    static constexpr long double skew = 12 * randlib::consts::sqrt2 * randlib::consts::sqrt3 * randlib::consts::apery / (randlib::consts::pi_sq * randlib::consts::pi);
  return skew;
}

template < typename RealType >
long double GumbelRand<RealType>::ExcessKurtosis() const
{
  return 2.4l;
}

template < typename RealType >
long double GumbelRand<RealType>::Entropy() const
{
  return logBeta + randlib::consts::egamma + 1.0;
}

} // namespace randlib

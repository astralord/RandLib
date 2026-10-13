#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The RaisedCosineDistribution class <BR>
 * Abstract class for Raised-cosine distribution
 *
 * Notation: X ~ Raised-cosine(μ, s)
 */
template <typename RealType>
class RaisedCosineDistribution : public ContinuousDistribution<RealType>
{
    double mu = 0; ///< location μ
    double s = randlib::consts::pi; ///< scale
    double s_pi = 1; ///< s / π
    double log2S = randlib::consts::ln2 + randlib::consts::lnpi; ///< log(2s)

protected:
    RaisedCosineDistribution(double location, double scale);

public:
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return mu - s; }
    [[nodiscard]] RealType MaxValue() const override { return mu + s; }

protected:
    void SetLocation(double location);
    void SetScale(double scale);

public:
    [[nodiscard]] inline double GetLocation() const noexcept { return mu; }
    [[nodiscard]] inline double GetScale() const noexcept { return s; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);
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
    std::complex<double> CFImpl(double t) const override;
};

/**
 * @brief The RaisedCosineRand class <BR>
 * Raised-cosine distribution
 */
template <typename RealType>
class RaisedCosineRand final : public RaisedCosineDistribution<RealType>
{
public:
    RaisedCosineRand(double location = 0, double scale = randlib::consts::pi) : RaisedCosineDistribution<RealType>(location, scale) {}
    [[nodiscard]] std::string Name() const override;

    using RaisedCosineDistribution<RealType>::SetLocation;
    using RaisedCosineDistribution<RealType>::SetScale;
};

/**
 * @brief The RaabGreenRand class <BR>
 * Raab-Green distribution
 *
 * Notation: X ~ Raab-Green()
 *
 * Related distributions:
 * X ~ Raised-cosine(0.0, π)
 */
template <typename RealType>
class RaabGreenRand final : public RaisedCosineDistribution<RealType>
{
public:
    RaabGreenRand() : RaisedCosineDistribution<RealType>(0.0, randlib::consts::pi) {}
    [[nodiscard]] std::string Name() const override;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template < typename RealType >
RaisedCosineDistribution<RealType>::RaisedCosineDistribution(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

template < typename RealType >
void RaisedCosineDistribution<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Raised cosine distribution: location should be finite");
    mu = location;
}

template < typename RealType >
void RaisedCosineDistribution<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Raised-Cosine distribution: scale should be positive");
    s = scale;
    s_pi = s * randlib::consts::inv_pi;
    log2S = std::log(2 * s);
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::f(const RealType &x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -static_cast<double>(randlib::consts::pi) || xAdj >= static_cast<double>(randlib::consts::pi))
        return 0.0;
    double y = std::cos(xAdj) + 1.0;
    return 0.5 * y / s;
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::logf(const RealType &x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -static_cast<double>(randlib::consts::pi) || xAdj >= static_cast<double>(randlib::consts::pi))
        return -INFINITY;
    double y = std::cos(xAdj);
    y = std::log1pl(y);
    return y - log2S;
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::F(const RealType &x) const
{
    double xAdj = (x - mu) / s;
    if (xAdj <= -1)
        return 0.0;
    if (xAdj >= 1)
        return 1.0;
    double y = std::sin(xAdj * randlib::consts::pi);
    y *= randlib::consts::inv_pi;
    y += xAdj + 1;
    return 0.5 * y;
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::S(const RealType &x) const
{
    double xAdj = (x - mu) / s;
    if (xAdj <= -1)
        return 1.0;
    if (xAdj >= 1)
        return 0.0;
    double y = std::sin(xAdj * randlib::consts::pi);
    y *= randlib::consts::inv_pi;
    y += xAdj;
    return 0.5 - 0.5 * y;
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType RaisedCosineDistribution<RealType>::StandardVariate(G &g)
{
    /// p. 160. Non-Uniform Random Variate Generation. Luc Devroye
    RealType X = randlib::consts::pi * detail::uniform01<RealType>(g) - randlib::consts::half_pi;
    RealType XSq = X * X;
    RealType U = 2 * detail::uniform01<RealType>(g);
    int a = 0, b = -1;
    RealType W = 0.0, V = 1.0;
    size_t iter = 0;
    do {
        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W += V;
        if (U >= W)
            return X;
        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W -= V;
        if (U <= W)
        {
            if (X == 0.0)
                return 0.0;
            return (X > 0.0) ? randlib::consts::pi - X : -randlib::consts::pi - X;
        }
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Raised-Cosine distribution: sampling failed");
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType RaisedCosineDistribution<RealType>::operator()(G &g) const
{
    return mu + s_pi * StandardVariate(g);
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::Mean() const
{
    return mu;
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::Variance() const
{
    static constexpr double coef = 1.0 / 3 - 2.0 / randlib::consts::pi_sq;
    return s * s * coef;
}

template < typename RealType >
std::complex<double> RaisedCosineDistribution<RealType>::CFImpl(double t) const
{
    /// π² sin(st) / (st (π² - (st)²)) = π² sinc(π - st) / (st (π + st)), without 0 / 0 at st = π
    double st = s * t;
    double u = static_cast<double>(randlib::consts::pi - st);
    double sinc = (u == 0.0) ? 1.0 : std::sin(u) / u;
    double y = randlib::consts::pi_sq * sinc / (st * (randlib::consts::pi + st));
    return y * std::exp(std::complex<double>(0.0, mu * t));
}

template < typename RealType >
RealType RaisedCosineDistribution<RealType>::Median() const
{
    return mu;
}

template < typename RealType >
RealType RaisedCosineDistribution<RealType>::Mode() const
{
    return mu;
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::Skewness() const
{
    return 0.0l;
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::ExcessKurtosis() const
{
    static constexpr long double numerator = 1.2 * (90.0 - randlib::consts::pi_sq * randlib::consts::pi_sq);
    static constexpr long double denominator = randlib::consts::pi_sq - 6.0;
    static constexpr long double y = numerator / (denominator * denominator);
    return y;
}
template < typename RealType >
std::string RaisedCosineRand<RealType>::Name() const
{
    return "Raised cosine(" + this->toStringWithPrecision(this->GetLocation()) + ", " + this->toStringWithPrecision(this->GetScale()) + ")";
}
template < typename RealType >
std::string RaabGreenRand<RealType>::Name() const
{
    return "Raab Green";
}

} // namespace randlib

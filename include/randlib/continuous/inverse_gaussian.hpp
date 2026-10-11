#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The InverseGaussianRand class <BR>
 * Inverse Gaussian (Wald) distribution
 *
 * Notation: X ~ IG(μ, λ)
 */
template <typename RealType>
class InverseGaussianRand final : public ContinuousDistribution<RealType>
{
    double mu = 1; ///< mean μ
    double lambda = 1; ///< shape λ
    double pdfCoef = randlib::consts::inv_sqrt2 / randlib::consts::sqrtpi; ///< (λ/(2π))^(1/2)
    double cdfCoef = randlib::consts::e * randlib::consts::e; ///< exp(2λ/μ)
public:
    InverseGaussianRand(double mean = 1, double shape = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double mean, double shape);
    [[nodiscard]] inline double GetMean() const noexcept { return mu; }
    [[nodiscard]] inline double GetShape() const noexcept { return lambda; }

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
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template < typename RealType >
InverseGaussianRand<RealType>::InverseGaussianRand(double mean, double shape)
{
    SetParameters(mean, shape);
}

template < typename RealType >
std::string InverseGaussianRand<RealType>::Name() const
{
    return "Inverse-Gaussian(" + this->toStringWithPrecision(GetMean()) + ", " + this->toStringWithPrecision(GetShape()) + ")";
}

template < typename RealType >
void InverseGaussianRand<RealType>::SetParameters(double mean, double shape)
{
    if (!(mean > 0.0))
        throw std::invalid_argument("Inverse-Gaussian distribution: mean should be positive");
    if (!(shape > 0.0))
        throw std::invalid_argument("Inverse-Gaussian distribution: shape should be positive");
    mu = mean;
    lambda = shape;

    pdfCoef = 0.5 * std::log(0.5 * lambda * randlib::consts::inv_pi);
    cdfCoef = std::exp(2 * lambda / mu);
}

template < typename RealType >
double InverseGaussianRand<RealType>::f(const RealType & x) const
{
    return (x > 0.0) ? std::exp(logf(x)) : 0.0;
}

template < typename RealType >
double InverseGaussianRand<RealType>::logf(const RealType & x) const
{
    if (x <= 0.0)
        return -INFINITY;
    double y = -1.5 * std::log(x);
    double z = (x - mu);
    z *= z;
    z *= -0.5 * lambda / (x * mu * mu);
    z += pdfCoef;
    return y + z;
}

template < typename RealType >
double InverseGaussianRand<RealType>::F(const RealType & x) const
{
    if (x <= 0.0)
        return 0.0;
    double b = std::sqrt(0.5 * lambda / x);
    double a = b * x / mu;
    if (b > a) {
        /// Lower tail. With 2λ/μ - (a + b)² = -(b - a)² both terms share the factor
        /// exp(-(b - a)²), so F = exp(-(b - a)²) (erfcx(b - a) + erfcx(a + b)) / 2 neither
        /// cancels nor overflows
        double d = b - a;
        return 0.5 * std::exp(-d * d) * (math::erfcx(d) + math::erfcx(a + b));
    }
    return 1.0 - S(x);
}

template < typename RealType >
double InverseGaussianRand<RealType>::S(const RealType & x) const
{
    if (x <= 0.0)
        return 1.0;
    double b = std::sqrt(0.5 * lambda / x);
    double a = b * x / mu;
    if (a >= b) {
        /// Upper tail: S = exp(-(a - b)²) (erfcx(a - b) - erfcx(a + b)) / 2, see F
        double d = a - b;
        return 0.5 * std::exp(-d * d) * (math::erfcx(d) - math::erfcx(a + b));
    }
    return 1.0 - F(x);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType InverseGaussianRand<RealType>::operator()(G &g) const
{
    RealType X = detail::standardNormal<RealType>(g);
    RealType U = detail::uniform01<RealType>(g);
    X *= X;
    RealType mupX = mu * X;
    RealType y = 4 * lambda + mupX;
    y = std::sqrt(y * mupX);
    y -= mupX;
    y *= -0.5 / lambda;
    ++y;
    if (U * (1 + y) > 1.0)
        y = 1.0 / y;
    return mu * y;
}

template < typename RealType >
long double InverseGaussianRand<RealType>::Mean() const
{
    return mu;
}

template < typename RealType >
long double InverseGaussianRand<RealType>::Variance() const
{
    return mu * mu * mu / lambda;
}

template < typename RealType >
std::complex<double> InverseGaussianRand<RealType>::CFImpl(double t) const
{
    double im = mu * mu;
    im *= t / lambda;
    std::complex<double> y(1, -im - im);
    y = 1.0 - std::sqrt(y);
    y *= lambda / mu;
    return std::exp(y);
}

template < typename RealType >
RealType InverseGaussianRand<RealType>::Mode() const
{
    RealType aux = 1.5 * mu / lambda;
    RealType mode = 1 + aux * aux;
    mode = std::sqrt(mode);
    mode -= aux;
    return mu * mode;
}

template < typename RealType >
long double InverseGaussianRand<RealType>::Skewness() const
{
    return 3 * std::sqrt(mu / lambda);
}

template < typename RealType >
long double InverseGaussianRand<RealType>::ExcessKurtosis() const
{
    return 15 * mu / lambda;
}

} // namespace randlib

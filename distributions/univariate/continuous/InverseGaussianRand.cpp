#include "InverseGaussianRand.h"
#include "NormalRand.h"
#include "UniformRand.h"

template < typename RealType >
InverseGaussianRand<RealType>::InverseGaussianRand(double mean, double shape)
{
    SetParameters(mean, shape);
}

template < typename RealType >
String InverseGaussianRand<RealType>::Name() const
{
    return "Inverse-Gaussian(" + this->toStringWithPrecision(GetMean()) + ", " + this->toStringWithPrecision(GetShape()) + ")";
}

template < typename RealType >
void InverseGaussianRand<RealType>::SetParameters(double mean, double shape)
{
    if (mean <= 0.0)
        throw std::invalid_argument("Inverse-Gaussian distribution: mean should be positive");
    if (shape <= 0.0)
        throw std::invalid_argument("Inverse-Gaussian distribution: shape should be positive");
    mu = mean;
    lambda = shape;

    pdfCoef = 0.5 * std::log(0.5 * lambda * M_1_PI);
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
    double y = std::erfc(b - a);
    y += cdfCoef * std::erfc(a + b);
    return 0.5 * y;
}

template < typename RealType >
double InverseGaussianRand<RealType>::S(const RealType & x) const
{
    if (x <= 0.0)
        return 1.0;
    double b = std::sqrt(0.5 * lambda / x);
    double a = b * x / mu;
    double y = std::erfc(a - b);
    y -= cdfCoef * std::erfc(a + b);
    return 0.5 * y;
}

template < typename RealType >
RealType InverseGaussianRand<RealType>::Variate() const
{
    RealType X = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType U = UniformRand::StandardVariate(this->localRandGenerator);
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

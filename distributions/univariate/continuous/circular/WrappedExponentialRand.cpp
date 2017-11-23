#include "WrappedExponentialRand.h"
#include "../UniformRand.h"

template < typename RealType >
WrappedExponentialRand<RealType>::WrappedExponentialRand(double rate) : CircularDistribution<RealType>(M_PI)
{
    SetRate(rate);
}

template < typename RealType >
String WrappedExponentialRand<RealType>::Name() const
{
    return "Wrapped Exponential(" + this->toStringWithPrecision(GetRate()) + ")";
}

template < typename RealType >
void WrappedExponentialRand<RealType>::SetRate(double rate)
{
    if (lambda <= 0.0)
        throw std::invalid_argument("Wrapped Exponential distribution: rate parameter should be positive");
    lambda = rate;
    logLambda = std::log(lambda);
    scaledLambda = 2 * M_PI * lambda;
    expmScaledLambda = std::exp(-scaledLambda);
    pdfCoef = -std::expm1l(-scaledLambda);
    logpdfCoef = RandMath::log1mexp(-scaledLambda);
}

template < typename RealType >
double WrappedExponentialRand<RealType>::f(const RealType &x) const
{
    return (x < 0 || x > 2 * M_PI) ? 0.0 : std::exp(logf(x));
}

template < typename RealType >
double WrappedExponentialRand<RealType>::logf(const RealType &x) const
{
    return (x < 0 || x > 2 * M_PI) ? - INFINITY : logLambda - lambda * x - logpdfCoef;
}

template < typename RealType >
double WrappedExponentialRand<RealType>::F(const RealType &x) const
{
    if (x <= 0.0)
        return 0.0;
    return (x < 2 * M_PI) ? std::exp(RandMath::log1mexp(-lambda * x) - logpdfCoef) : 1.0;
}

template < typename RealType >
double WrappedExponentialRand<RealType>::S(const RealType &x) const
{
    if (x <= 0.0)
        return 1.0;
    if (x >= 2 * M_PI)
        return 0.0;
    double y = std::expm1l(scaledLambda - lambda * x);
    y /= pdfCoef;
    return expmScaledLambda * y;
}

template < typename RealType >
RealType WrappedExponentialRand<RealType>::Variate() const
{
    return quantileImpl(UniformRand<RealType>::StandardVariate(this->localRandGenerator));
}

template < typename RealType >
long double WrappedExponentialRand<RealType>::CircularMean() const
{
    return M_PI_2 - RandMath::atan(lambda);
}

template < typename RealType >
long double WrappedExponentialRand<RealType>::CircularVariance() const
{
    return 1.0 - 1.0 / std::sqrt(1.0 + lambda * lambda);
}

template < typename RealType >
RealType WrappedExponentialRand<RealType>::Median() const
{
    return (M_LN2 - RandMath::log1pexp(-scaledLambda)) / lambda;
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
    double temp = t / lambda;
    double coef = 1.0 / (1.0 + temp * temp);
    return std::complex<double>(coef, temp * coef);
}


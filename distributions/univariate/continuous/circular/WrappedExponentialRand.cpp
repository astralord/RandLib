#include "WrappedExponentialRand.h"
#include "../UniformRand.h"

template < typename RealType >
WrappedExponentialRand<RealType>::WrappedExponentialRand(double rate) : CircularDistribution<RealType>(randlib::consts::pi)
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
    if (!(rate > 0.0))
        throw std::invalid_argument("Wrapped Exponential distribution: rate parameter should be positive");
    lambda = rate;
    logLambda = std::log(lambda);
    scaledLambda = 2 * randlib::consts::pi * lambda;
    expmScaledLambda = std::exp(-scaledLambda);
    pdfCoef = -std::expm1l(-scaledLambda);
    logpdfCoef = RandMath::log1mexp(-scaledLambda);
}

template < typename RealType >
double WrappedExponentialRand<RealType>::f(const RealType &x) const
{
    return (x < this->MinValue() || x > this->MaxValue()) ? 0.0 : std::exp(logf(x));
}

template < typename RealType >
double WrappedExponentialRand<RealType>::logf(const RealType &x) const
{
    return (x < this->MinValue() || x > this->MaxValue()) ? - INFINITY : logLambda - lambda * x - logpdfCoef;
}

template < typename RealType >
double WrappedExponentialRand<RealType>::F(const RealType &x) const
{
    if (x <= 0.0)
        return 0.0;
    return (x < this->MaxValue()) ? std::exp(RandMath::log1mexp(-lambda * x) - logpdfCoef) : 1.0;
}

template < typename RealType >
double WrappedExponentialRand<RealType>::S(const RealType &x) const
{
    if (x <= 0.0)
        return 1.0;
    if (x >= this->MaxValue())
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
    return randlib::consts::half_pi - RandMath::atan(lambda);
}

template < typename RealType >
long double WrappedExponentialRand<RealType>::CircularVariance() const
{
    return 1.0 - lambda / std::sqrt(1.0 + lambda * lambda);
}

template < typename RealType >
RealType WrappedExponentialRand<RealType>::Median() const
{
    return (randlib::consts::ln2 - RandMath::softplus(-scaledLambda)) / lambda;
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
    /// φ(t) = λ (1 - e^{2π(it - λ)}) / ((λ - it)(1 - e^{-2πλ})), which reduces to λ / (λ - it) for integer t
    std::complex<double> lambdaMinusIt(lambda, -t);
    std::complex<double> wrap = 1.0 - std::exp(-lambdaMinusIt * static_cast<double>(2 * randlib::consts::pi));
    return lambda * wrap / (lambdaMinusIt * pdfCoef);
}


template class WrappedExponentialRand<float>;
template class WrappedExponentialRand<double>;
template class WrappedExponentialRand<long double>;

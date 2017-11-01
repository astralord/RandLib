#include "WrappedExponentialRand.h"
#include "../UniformRand.h"

WrappedExponentialRand::WrappedExponentialRand(double rate) : CircularDistribution(M_PI)
{
    SetRate(rate);
}

String WrappedExponentialRand::Name() const
{
    return "Wrapped Exponential(" + toStringWithPrecision(GetRate()) + ")";
}

void WrappedExponentialRand::SetRate(double rate)
{
    if (lambda <= 0.0)
        throw std::invalid_argument("Wrapped Exponential distribution: rate parameter should be positive");
    lambda = rate;
    logLambda = std::log(lambda);
    scaledLambda = 2 * M_PI * lambda;
    expmScaledLambda = std::exp(-scaledLambda);
    pdfCoef = -std::expm1(-scaledLambda);
    logpdfCoef = RandMath::log1mexp(-scaledLambda);
}

double WrappedExponentialRand::f(const double &x) const
{
    return (x < 0 || x > 2 * M_PI) ? 0.0 : std::exp(logf(x));
}

double WrappedExponentialRand::logf(const double &x) const
{
    return (x < 0 || x > 2 * M_PI) ? - INFINITY : logLambda - lambda * x - logpdfCoef;
}

double WrappedExponentialRand::F(const double &x) const
{
    if (x <= 0.0)
        return 0.0;
    return (x < 2 * M_PI) ? std::exp(RandMath::log1mexp(-lambda * x) - logpdfCoef) : 1.0;
}

double WrappedExponentialRand::S(const double &x) const
{
    if (x <= 0.0)
        return 1.0;
    if (x >= 2 * M_PI)
        return 0.0;
    double y = std::expm1(scaledLambda - lambda * x);
    y /= pdfCoef;
    return expmScaledLambda * y;
}

double WrappedExponentialRand::Variate() const
{
    return quantileImpl(UniformRand::StandardVariate());
}

double WrappedExponentialRand::CircularMean() const
{
    return M_PI_2 - RandMath::atan(lambda);
}

double WrappedExponentialRand::CircularVariance() const
{
    return 1.0 - 1.0 / std::sqrt(1.0 + lambda * lambda);
}

double WrappedExponentialRand::Median() const
{
    return (M_LN2 - RandMath::log1pexp(-scaledLambda)) / lambda;
}

double WrappedExponentialRand::Mode() const
{
    return 0.0;
}

double WrappedExponentialRand::quantileImpl(double p) const
{
    return -std::log1p(-p * pdfCoef) / lambda;
}

double WrappedExponentialRand::quantileImpl1m(double p) const
{
    return -std::log(expmScaledLambda + p * pdfCoef) / lambda;
}

std::complex<double> WrappedExponentialRand::CFImpl(double t) const
{
    double temp = t / lambda;
    double coef = 1.0 / (1.0 + temp * temp);
    return std::complex<double>(coef, temp * coef);
}


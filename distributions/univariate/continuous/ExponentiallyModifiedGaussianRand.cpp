#include "ExponentiallyModifiedGaussianRand.h"

template < typename RealType >
ExponentiallyModifiedGaussianRand<RealType>::ExponentiallyModifiedGaussianRand(double location, double variance, double rate)
{
    SetParameters(location, variance, rate);
}

template < typename RealType >
String ExponentiallyModifiedGaussianRand<RealType>::Name() const
{
    return "Exponentially modified Gaussian(" + this->toStringWithPrecision(GetLocation()) + ", "
                                              + this->toStringWithPrecision(X.Variance()) + ", "
                                              + this->toStringWithPrecision(GetRate()) + ")";
}

template < typename RealType >
void ExponentiallyModifiedGaussianRand<RealType>::SetParameters(double location, double variance, double rate)
{
    if (variance <= 0)
        throw std::invalid_argument("Exponentially modified Gaussian distribution: variance should be positive");
    if (rate <= 0)
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
    b = M_SQRT1_2 / sigma;
    v = beta * sigma;
}

template < typename RealType >
double ExponentiallyModifiedGaussianRand<RealType>::f(const RealType & x) const
{
    double lambda = Y.GetRate();
    double y = a - x;
    y *= b;
    y = std::erfc(y);
    y *= 0.5 * lambda;
    double exponent = c - x;
    exponent *= lambda;
    exponent = std::exp(exponent);
    return y * exponent;
}

template < typename RealType >
double ExponentiallyModifiedGaussianRand<RealType>::logf(const RealType & x) const
{
    double lambda = Y.GetRate();
    double y = a - x;
    y *= b;
    y = std::erfc(y);
    y *= 0.5 * lambda;
    double exponent = c - x;
    exponent *= lambda;
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
RealType ExponentiallyModifiedGaussianRand<RealType>::Variate() const
{
    return X.Variate() + Y.Variate();
}

template < typename RealType >
RealType ExponentiallyModifiedGaussianRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    return NormalRand<RealType>::StandardVariate(randGenerator) + ExponentialRand::StandardVariate(randGenerator);
}

template < typename RealType >
void ExponentiallyModifiedGaussianRand<RealType>::Reseed(unsigned long seed) const
{
    X.Reseed(seed);
    Y.Reseed(seed + 1);
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


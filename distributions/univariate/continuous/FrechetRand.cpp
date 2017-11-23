#include "FrechetRand.h"
#include "ExponentialRand.h"

template < typename RealType >
FrechetRand<RealType>::FrechetRand(double shape, double scale, double location)
{
    SetParameters(shape, scale, location);
}

template < typename RealType >
String FrechetRand<RealType>::Name() const
{
    return "Frechet(" + this->toStringWithPrecision(GetShape()) + ", "
                      + this->toStringWithPrecision(GetScale()) + ", "
                      + this->toStringWithPrecision(GetLocation()) + ")";
}

template < typename RealType >
void FrechetRand<RealType>::SetParameters(double shape, double scale, double location)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Frechet distribution: shape should be positive");
    if (scale <= 0.0)
        throw std::invalid_argument("Frechet distribution: scale should be positive");
    alpha = shape;
    alphaInv = 1.0 / alpha;
    s = scale;
    pdfCoef = std::log(alpha / s);
    m = location;
}

template < typename RealType >
double FrechetRand<RealType>::f(const RealType &x) const
{
    return (x <= m) ? 0.0 : std::exp(logf(x));
}

template < typename RealType >
double FrechetRand<RealType>::logf(const RealType &x) const
{
    if (x <= m)
        return -INFINITY;
    double xAdj = (x - m) / s;
    double logxAdj = std::log(xAdj);
    double a = alpha * logxAdj;
    double expA = std::exp(-a);
    return pdfCoef - a - expA - logxAdj;
}

template < typename RealType >
double FrechetRand<RealType>::F(const RealType & x) const
{
    if (x <= m)
        return 0.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
    return std::exp(-xPow);
}

template < typename RealType >
double FrechetRand<RealType>::S(const RealType & x) const
{
    if (x <= m)
        return 1.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
    return -std::expm1l(-xPow);
}

template < typename RealType >
RealType FrechetRand<RealType>::Variate() const
{
    return m + s / std::pow(ExponentialRand<RealType>::StandardVariate(this->localRandGenerator), alphaInv);
}

template < typename RealType >
long double FrechetRand<RealType>::Mean() const
{
    if (alpha <= 1.0)
        return INFINITY;
    return m + s * std::tgammal(1.0 - alphaInv);
}

template < typename RealType >
long double FrechetRand<RealType>::Variance() const
{
    if (alpha <= 2.0)
        return INFINITY;
    long double var = std::tgammal(1.0 - alphaInv);
    var *= var;
    var = std::tgammal(1.0 - 2 * alphaInv) - var;
    return s * s * var;
}

template < typename RealType >
RealType FrechetRand<RealType>::quantileImpl(double p) const
{
    RealType y = -std::log(p);
    y = s / std::pow(y, alphaInv);
    return y + m;
}

template < typename RealType >
RealType FrechetRand<RealType>::quantileImpl1m(double p) const
{
    RealType y = -std::log1pl(-p);
    y = s / std::pow(y, alphaInv);
    return y + m;
}

template < typename RealType >
RealType FrechetRand<RealType>::Median() const
{
    return m + s / std::pow(M_LN2, alphaInv);
}

template < typename RealType >
RealType FrechetRand<RealType>::Mode() const
{
    RealType y = alpha / (1.0 + alpha);
    y = std::pow(y, alphaInv);
    return m + s * y;
}

template < typename RealType >
long double FrechetRand<RealType>::Skewness() const
{
    if (alpha <= 3.0)
        return INFINITY;
    long double x = std::tgammal(1.0 - alphaInv);
    long double y = std::tgammal(1.0 - 2.0 * alphaInv);
    long double z = std::tgammal(1.0 - 3.0 * alphaInv);
    long double numerator = 2 * x * x - 3 * y;
    numerator *= x;
    numerator += z;
    long double denominator = y - x * x;
    denominator = std::pow(denominator, 1.5);
    return numerator / denominator;
}

template < typename RealType >
long double FrechetRand<RealType>::ExcessKurtosis() const
{
    if (alpha <= 4.0)
        return INFINITY;
    long double x = std::tgammal(1.0 - alphaInv);
    long double y = std::tgammal(1.0 - 2.0 * alphaInv);
    long double z = std::tgammal(1.0 - 3.0 * alphaInv);
    long double w = std::tgammal(1.0 - 4.0 * alphaInv);
    long double numerator = w - 4 * z * x + 3 * y * y;
    long double denominator = y - x * x;
    denominator *= denominator;
    return numerator / denominator - 6.0;
}

template < typename RealType >
double FrechetRand<RealType>::Entropy() const
{
    return 1.0 + M_EULER * (1.0 + alphaInv) + std::log(s / alpha);
}

template class FrechetRand<float>;
template class FrechetRand<double>;
template class FrechetRand<long double>;

#include "YuleRand.h"

template < typename IntType >
YuleRand<IntType>::YuleRand(double shape) :
X(shape, 1.0)
{
    SetShape(shape);
}

template < typename IntType >
String YuleRand<IntType>::Name() const
{
    return "Yule(" + this->toStringWithPrecision(GetShape()) + ")";
}

template < typename IntType >
void YuleRand<IntType>::SetShape(double shape)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Yule distribution: shape should be positive");
    rho = shape;
    lgamma1pRo = std::lgammal(rho + 1);
    X.SetShape(rho);
}

template < typename IntType >
double YuleRand<IntType>::logP(const IntType & k) const
{
    if (k < 1)
        return -INFINITY;
    double y = lgamma1pRo;
    y += RandMath::lfact(k - 1);
    y -= std::lgammal(k + rho + 1);
    y += X.GetLogShape();
    return y;
}

template < typename IntType >
double YuleRand<IntType>::F(const IntType & k) const
{
    if (k < 1)
        return 0.0;
    double y = lgamma1pRo;
    y += RandMath::lfact(k - 1);
    y -= std::lgammal(k + rho + 1);
    double logk = std::log(k);
    return -std::expm1(y + logk);
}

template < typename IntType >
double YuleRand<IntType>::S(const IntType & k) const
{
    if (k < 1)
        return 1.0;
    double y = lgamma1pRo;
    y += RandMath::lfact(k - 1);
    y -= std::lgammal(k + rho + 1);
    y = std::exp(y);
    return k * y;
}

template < typename IntType >
IntType YuleRand<IntType>::Variate() const
{
    double prob = 1.0 / X.Variate();
    return GeometricRand<IntType>::Variate(prob, this->localRandGenerator) + 1;
}

template < typename IntType >
IntType YuleRand<IntType>::Variate(double shape, RandGenerator &randGenerator)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Yule distribution: shape should be positive");
    double prob = 1.0 / ParetoRand<double>::StandardVariate(shape, randGenerator);
    return GeometricRand<IntType>::Variate(prob, randGenerator) + 1;
}

template < typename IntType >
void YuleRand<IntType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    X.Reseed(seed + 1);
}

template < typename IntType >
long double YuleRand<IntType>::Mean() const
{
    return (rho <= 1) ? INFINITY : rho / (rho - 1);
}

template < typename IntType >
long double YuleRand<IntType>::Variance() const
{
    if (rho <= 2)
        return INFINITY;
    double aux = rho / (rho - 1);
    return aux * aux / (rho - 2);
}

template < typename IntType >
IntType YuleRand<IntType>::Mode() const
{
    return 1;
}

template < typename IntType >
long double YuleRand<IntType>::Skewness() const
{
    if (rho <= 3)
        return INFINITY;
    long double skewness = rho + 1;
    skewness *= skewness;
    skewness *= std::sqrt(rho - 2);
    return skewness / (rho * (rho - 3));
}

template < typename IntType >
long double YuleRand<IntType>::ExcessKurtosis() const
{
    if (rho <= 4)
        return INFINITY;
    long double numerator = 11 * rho * rho - 49;
    numerator *= rho;
    numerator -= 22;
    long double denominator = rho * (rho - 4) * (rho - 3);
    return rho + 3 + numerator / denominator;
}


template class YuleRand<int>;
template class YuleRand<long int>;
template class YuleRand<long long int>;

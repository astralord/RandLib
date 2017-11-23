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
    ro = shape;
    lgamma1pRo = std::lgammal(ro + 1);
    X.SetShape(ro);
}

template < typename IntType >
double YuleRand<IntType>::P(const IntType & k) const
{
    return (k < 1) ? 0.0 : std::exp(logP(k));
}

template < typename IntType >
double YuleRand<IntType>::logP(const IntType & k) const
{
    if (k < 1)
        return -INFINITY;
    double y = lgamma1pRo;
    y += RandMath::lfact(k - 1);
    y -= std::lgammal(k + ro + 1);
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
    y -= std::lgammal(k + ro + 1);
    y = std::exp(y);
    return 1.0 - k * y;
}

template < typename IntType >
double YuleRand<IntType>::S(const IntType & k) const
{
    if (k < 1)
        return 1.0;
    double y = lgamma1pRo;
    y += RandMath::lfact(k - 1);
    y -= std::lgammal(k + ro + 1);
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
        return -1;
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
    return (ro <= 1) ? INFINITY : ro / (ro - 1);
}

template < typename IntType >
long double YuleRand<IntType>::Variance() const
{
    if (ro <= 2)
        return INFINITY;
    double aux = ro / (ro - 1);
    return aux * aux / (ro - 2);
}

template < typename IntType >
IntType YuleRand<IntType>::Mode() const
{
    return 1;
}

template < typename IntType >
long double YuleRand<IntType>::Skewness() const
{
    if (ro <= 3)
        return INFINITY;
    long double skewness = ro + 1;
    skewness *= skewness;
    skewness *= std::sqrt(ro - 2);
    return skewness / (ro * (ro - 3));
}

template < typename IntType >
long double YuleRand<IntType>::ExcessKurtosis() const
{
    if (ro <= 4)
        return INFINITY;
    long double numerator = 11 * ro * ro - 49;
    numerator *= ro;
    numerator -= 22;
    long double denominator = ro * (ro - 4) * (ro - 3);
    return ro + 3 + numerator / denominator;
}


template class YuleRand<int>;
template class YuleRand<long int>;
template class YuleRand<long long int>;

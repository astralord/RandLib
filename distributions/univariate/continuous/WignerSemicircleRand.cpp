#include "WignerSemicircleRand.h"

template < typename RealType >
WignerSemicircleRand<RealType>::WignerSemicircleRand(double radius)
{
    SetRadius(radius);
}

template < typename RealType >
String WignerSemicircleRand<RealType>::Name() const
{
    return "Wigner Semicircle(" + this->toStringWithPrecision(GetRadius()) + ")";
}

template < typename RealType >
void WignerSemicircleRand<RealType>::SetRadius(double radius)
{
    if (radius <= 0.0)
        throw std::invalid_argument("Wigner-Semicircle distribution: radius should be positive");
    R = radius;
    RSq = R * R;
    logRSq = std::log(RSq);
}

template < typename RealType >
double WignerSemicircleRand<RealType>::f(const RealType & x) const
{
    double xSq = x * x;
    if (xSq >= RSq)
        return 0.0;
    double y = RSq - xSq;
    y = std::sqrt(y);
    y *= M_1_PI / RSq;
    return 2 * y;
}

template < typename RealType >
double WignerSemicircleRand<RealType>::logf(const RealType & x) const
{
    double xSq = x * x;
    if (xSq >= RSq)
        return -INFINITY;
    return M_LN2 + 0.5 * std::log(RSq - xSq) - M_LNPI - logRSq;
}

template < typename RealType >
double WignerSemicircleRand<RealType>::F(const RealType &x) const
{
    if (x <= -R)
        return 0.0;
    if (x >= R)
        return 1.0;
    double y = RSq - x * x;
    y = x * std::sqrt(y) / RSq;
    double z = std::asin(x / R);
    return 0.5 + (y + z) / M_PI;
}

template < typename RealType >
RealType WignerSemicircleRand<RealType>::Variate() const
{
    RealType x = X.Variate();
    x += x - 1;
    return R * x;
}

template < typename RealType >
void WignerSemicircleRand<RealType>::Reseed(unsigned long seed) const
{
    X.Reseed(seed);
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::Mean() const
{
    return 0.0;
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::Variance() const
{
    return 0.25 * RSq;
}

template < typename RealType >
RealType WignerSemicircleRand<RealType>::Median() const
{
    return 0.0;
}

template < typename RealType >
RealType WignerSemicircleRand<RealType>::Mode() const
{
    return 0.0;
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::Skewness() const
{
    return 0.0;
}

template < typename RealType >
long double WignerSemicircleRand<RealType>::ExcessKurtosis() const
{
    return -1.0;
}

template < typename RealType >
double WignerSemicircleRand<RealType>::Entropy() const
{
    return M_LNPI + 0.5 * logRSq - 0.5;
}

template class WignerSemicircleRand<float>;
template class WignerSemicircleRand<double>;
template class WignerSemicircleRand<long double>;

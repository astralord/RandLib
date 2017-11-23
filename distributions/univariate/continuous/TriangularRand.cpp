#include "TriangularRand.h"
#include "UniformRand.h"

template < typename RealType >
TriangularRand<RealType>::TriangularRand(double lowerLimit, double mode, double upperLimit)
{
    SetParameters(lowerLimit, mode, upperLimit);
}

template < typename RealType >
String TriangularRand<RealType>::Name() const
{
    return "Triangular("
            + this->toStringWithPrecision(MinValue()) + ", "
            + this->toStringWithPrecision(Mode()) + ", "
            + this->toStringWithPrecision(MaxValue()) + ")";
}

template < typename RealType >
void TriangularRand<RealType>::SetParameters(double lowerLimit, double mode, double upperLimit)
{
    if (lowerLimit >= mode)
        throw std::invalid_argument("Triangular distribution: lower limit should be larger than mode");
    if (mode >= upperLimit)
        throw std::invalid_argument("Triangular distribution: upper limit should be smaller than mode");
    a = lowerLimit;
    c = mode;
    b = upperLimit;
    SetConstantsForGenerator();
}

template < typename RealType >
void TriangularRand<RealType>::SetConstantsForGenerator()
{
    constForGenerator = (c - a) / (b - a);
    coefGenerator1 = (b - a) * (c - a);
    coefGenerator2 = (b - a) * (b - c);
}

template < typename RealType >
double TriangularRand<RealType>::f(const RealType & x) const
{
    if (x <= a)
        return 0;
    if (x < c)
        return 2.0 * (x - a) / ((b - a) * (c - a));
    if (x == c)
        return 2.0 / (b - a);
    if (x < b)
        return 2.0 * (b - x) / ((b - a) * (b - c));
    return 0;
}

template < typename RealType >
double TriangularRand<RealType>::logf(const RealType & x) const
{
    return std::log(f(x));
}

template < typename RealType >
double TriangularRand<RealType>::F(const RealType & x) const
{
    if (x <= a)
        return 0.0;
    if (x <= c)
        return (x - a) * (x - a) / ((b - a) * (c - a));
    if (x < b)
        return 1.0 - (b - x) * (b - x) / ((b - a) * (b - c));
    return 1.0;
}

template < typename RealType >
double TriangularRand<RealType>::S(const RealType & x) const
{
    if (x <= a)
        return 1.0;
    if (x <= c)
        return 1.0 - (x - a) * (x - a) / ((b - a) * (c - a));
    if (x < b)
        return (b - x) * (b - x) / ((b - a) * (b - c));
    return 0.0;
}

template < typename RealType >
RealType TriangularRand<RealType>::Variate() const
{
    RealType U = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
    if (U < constForGenerator)
        return a + std::sqrt(U * coefGenerator1);
    return b - std::sqrt((1 - U) * coefGenerator2);
}

template < typename RealType >
long double TriangularRand<RealType>::Mean() const
{
    return (a + b + c) / 3.0;
}

template < typename RealType >
long double TriangularRand<RealType>::Variance() const
{
    return (a * (a - b) + b * (b - c) + c * (c - a)) / 18.0;
}

template < typename RealType >
std::complex<double> TriangularRand<RealType>::CFImpl(double t) const
{
    double bmc = b - c, bma = b - a, cma = c - a;
    double at = a * t, bt = b * t, ct = c * t;
    std::complex<double> x(bmc * std::cos(at), bmc * std::sin(at));
    std::complex<double> y(bma * std::cos(ct), bma * std::sin(ct));
    std::complex<double> z(cma * std::cos(bt), cma * std::sin(bt));
    std::complex<double> numerator = x - y + z;
    /// in order to avoid numerical errors
    if (t < 1e-10 && std::fabs(numerator.real()) < 1e-10)
        return 1;
    double denominator = bma * cma * bmc * t * t;
    std::complex<double> frac = -numerator / denominator;
    return frac + frac;
}

template < typename RealType >
RealType TriangularRand<RealType>::Median() const
{
    if (c + c > a + b)
        return a + std::sqrt(0.5 * (b - a) * (c - a));
    return b - std::sqrt(0.5 * (b - a) * (b - c));
}

template < typename RealType >
RealType TriangularRand<RealType>::Mode() const
{
    return c;
}

template < typename RealType >
long double TriangularRand<RealType>::Skewness() const
{
    double numerator = M_SQRT2;
    numerator *= (a + b - c - c);
    numerator *= (a + a - b - c);
    numerator *= (a - b - b + c);
    double denominator = a * (a - b);
    denominator += b * (b - c);
    denominator += c * (c - a);
    denominator *= std::sqrt(denominator);
    return 0.2 * numerator / denominator;
}

template < typename RealType >
long double TriangularRand<RealType>::ExcessKurtosis() const
{
    return -0.6;
}


template class TriangularRand<float>;
template class TriangularRand<double>;
template class TriangularRand<long double>;

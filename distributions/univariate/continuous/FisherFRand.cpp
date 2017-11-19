#include "FisherFRand.h"

template < typename RealType >
FisherFRand<RealType>::FisherFRand(int degree1, int degree2)
{
    SetDegrees(degree1, degree2);
}

template < typename RealType >
String FisherFRand<RealType>::Name() const
{
    return "Fisher-F(" + this->toStringWithPrecision(GetFirstDegree()) + ", " + this->toStringWithPrecision(GetSecondDegree()) + ")";
}

template < typename RealType >
void FisherFRand<RealType>::SetDegrees(int degree1, int degree2)
{
    if (degree1 <= 0 || degree2 <= 0)
        throw std::invalid_argument("F-distribution: degrees of should be positive");

    d1 = degree1;
    d2 = degree2;

    B.SetShapes(0.5 * d1, 0.5 * d2);

    a = 0.5 * d1 - 1;
    d1_d2 = static_cast<double>(d1) / d2;
    c = -0.5 * (d1 + d2);
    d2_d1 = 1.0 / d1_d2;

    pdfCoef = (a + 1) * std::log(d1_d2);
    pdfCoef -= B.GetLogBetaFunction();
}

template < typename RealType >
double FisherFRand<RealType>::f(const RealType & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (a == 0.0)
            return std::exp(pdfCoef);
        return (a > 0) ? 0.0 : INFINITY;
    }
    return std::exp(logf(x));
}

template < typename RealType >
double FisherFRand<RealType>::logf(const RealType & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0) {
        if (a == 0.0)
            return pdfCoef;
        return (a > 0) ? -INFINITY : INFINITY;
    }
    double y = a * std::log(x);
    y += c * std::log1pl(d1_d2 * x);
    return pdfCoef + y;
}

template < typename RealType >
double FisherFRand<RealType>::F(const RealType & x) const
{
    return B.F(d1_d2 * x);
}

template < typename RealType >
double FisherFRand<RealType>::S(const RealType & x) const
{
    return B.S(d1_d2 * x);
}

template < typename RealType >
RealType FisherFRand<RealType>::Variate() const
{
    return d2_d1 * B.Variate();
}

template < typename RealType >
void FisherFRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    B.Sample(outputData);
    for (RealType &var : outputData)
        var = d2_d1 * var;
}

template < typename RealType >
void FisherFRand<RealType>::Reseed(unsigned long seed) const
{
    B.Reseed(seed);
}

template < typename RealType >
long double FisherFRand<RealType>::Mean() const
{
    return (d2 > 2) ? 1 + 2.0 / (d2 - 2) : INFINITY;
}

template < typename RealType >
long double FisherFRand<RealType>::Variance() const
{
    if (d2 <= 4)
        return INFINITY;
    double variance = d2;
    variance /= d2 - 2;
    variance *= variance;
    variance *= 2 * (d1 + d2 - 2);
    variance /= d1;
    variance /= d2 - 4;
    return variance;
}

template < typename RealType >
RealType FisherFRand<RealType>::Median() const
{
    return d2_d1 * B.Median();
}

template < typename RealType >
RealType FisherFRand<RealType>::Mode() const
{
    if (d1 <= 2)
        return 0.0;
    return d2_d1 * (d1 - 2) / (d2 + 2);
}

template < typename RealType >
long double FisherFRand<RealType>::Skewness() const
{
    if (d2 <= 6)
        return INFINITY;
    long double skewness = 8.0 * (d2 - 4.0);
    long double aux = d1 + d2 - 2;
    skewness /= d1 * aux;
    skewness = std::sqrt(skewness);
    skewness *= d1 + aux;
    skewness /= d2 - 6.0;
    return skewness;
}

template < typename RealType >
long double FisherFRand<RealType>::ExcessKurtosis() const
{
    if (d2 <= 8)
        return INFINITY;
    long double kurtosis = d2 - 2;
    kurtosis *= kurtosis;
    kurtosis *= d2 - 4;
    kurtosis /= d1;
    kurtosis /= d1 + d2 - 2;
    kurtosis += 5 * d2 - 22;
    kurtosis /= d2 - 6;
    kurtosis /= d2 - 8;
    return 12.0 * kurtosis;
}

template < typename RealType >
RealType FisherFRand<RealType>::quantileImpl(double p) const
{
  return d2_d1 * B.Quantile(p);
}

template < typename RealType >
RealType FisherFRand<RealType>::quantileImpl1m(double p) const
{
  return d2_d1 * B.Quantile1m(p);
}

template < typename RealType >
std::complex<double> FisherFRand<RealType>::CFImpl(double t) const
{
    return B.CF(d2_d1 * t);
}

template class FisherFRand<float>;
template class FisherFRand<double>;
template class FisherFRand<long double>;

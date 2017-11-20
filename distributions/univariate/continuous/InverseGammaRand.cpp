#include "InverseGammaRand.h"

template < typename RealType >
InverseGammaRand<RealType>::InverseGammaRand(double shape, double rate)
{
    SetParameters(shape, rate);
}

template < typename RealType >
String InverseGammaRand<RealType>::Name() const
{
    return "Inverse-Gamma(" + this->toStringWithPrecision(GetShape()) + ", " + this->toStringWithPrecision(GetRate()) + ")";
}

template < typename RealType >
void InverseGammaRand<RealType>::SetParameters(double shape, double rate)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Inverse-Gamma distribution: shape should be positive");
    if (rate <= 0.0)
        throw std::invalid_argument("Inverse-Gamma distribution: rate should be positive");
    X.SetParameters(shape, rate);
    alpha = X.GetShape();
    beta = X.GetRate();
    pdfCoef = -X.GetLogGammaShape() + alpha * X.GetLogRate();
}

template < typename RealType >
double InverseGammaRand<RealType>::f(const RealType & x) const
{
    return (x > 0.0) ? std::exp(logf(x)) : 0.0;
}

template < typename RealType >
double InverseGammaRand<RealType>::logf(const RealType & x) const
{
    if (x <= 0.0)
        return -INFINITY;
    double logX = std::log(x);
    double y = -(alpha - 1.0) * logX;
    y -= beta / x;
    y += pdfCoef;
    return y - 2 * logX;
}

template < typename RealType >
double InverseGammaRand<RealType>::F(const RealType & x) const
{
    return (x > 0.0) ? X.S(1.0 / x) : 0.0;
}

template < typename RealType >
double InverseGammaRand<RealType>::S(const RealType & x) const
{
    return (x > 0.0) ? X.F(1.0 / x) : 1.0;
}

template < typename RealType >
RealType InverseGammaRand<RealType>::Variate() const
{
    return 1.0 / X.Variate();
}

template < typename RealType >
void InverseGammaRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    X.Sample(outputData);
    for (RealType &var : outputData)
        var = 1.0 / var;
}

template < typename RealType >
void InverseGammaRand<RealType>::Reseed(unsigned long seed) const
{
    X.Reseed(seed);
}

template < typename RealType >
long double InverseGammaRand<RealType>::Mean() const
{
    return (alpha > 1) ? beta / (alpha - 1) : INFINITY;
}

template < typename RealType >
long double InverseGammaRand<RealType>::Variance() const
{
    if (alpha <= 2)
        return INFINITY;
    double var = beta / (alpha - 1);
    var *= var;
    return var / (alpha - 2);
}

template < typename RealType >
RealType InverseGammaRand<RealType>::Median() const
{
    return 1.0 / X.Median();
}

template < typename RealType >
RealType InverseGammaRand<RealType>::Mode() const
{
    return beta / (alpha + 1);
}


template < typename RealType >
RealType InverseGammaRand<RealType>::quantileImpl(double p) const
{
    return 1.0 / X.Quantile1m(p);
}

template < typename RealType >
RealType InverseGammaRand<RealType>::quantileImpl1m(double p) const
{
    return 1.0 / X.Quantile(p);
}

template < typename RealType >
long double InverseGammaRand<RealType>::Skewness() const
{
    return (alpha > 3) ? 4 * std::sqrt(alpha - 2) / (alpha - 3) : INFINITY;
}

template < typename RealType >
long double InverseGammaRand<RealType>::ExcessKurtosis() const
{
    if (alpha <= 4)
        return INFINITY;
    long double numerator = 30 * alpha - 66.0;
    long double denominator = (alpha - 3) * (alpha - 4);
    return numerator / denominator;
}

template class InverseGammaRand<float>;
template class InverseGammaRand<double>;
template class InverseGammaRand<long double>;

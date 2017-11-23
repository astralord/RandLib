#include "SkellamRand.h"
#include "../continuous/NoncentralChiSquaredRand.h"

template < typename IntType >
SkellamRand<IntType>::SkellamRand(double rate1, double rate2)
{
    SetRates(rate1, rate2);
}

template < typename IntType >
String SkellamRand<IntType>::Name() const
{
    return "Skellam(" + this->toStringWithPrecision(GetFirstRate()) + ", " + this->toStringWithPrecision(GetSecondRate()) + ")";
}

template < typename IntType >
void SkellamRand<IntType>::SetRates(double rate1, double rate2)
{
    if (rate1 <= 0.0 || rate2 <= 0.0)
        throw std::invalid_argument("Skellam distribution: rates should be positive");

    X.SetRate(rate1);
    mu1 = X.GetRate();
    logMu1 = std::log(mu1);
    sqrtMu1 = std::sqrt(mu1);

    Y.SetRate(rate2);
    mu2 = Y.GetRate();
    logMu2 = std::log(mu2);
    sqrtMu2 = std::sqrt(mu2);
}

template < typename IntType >
double SkellamRand<IntType>::P(const IntType & k) const
{
    return std::exp(logP(k));
}

template < typename IntType >
double SkellamRand<IntType>::logP(const IntType & k) const
{
    double y = RandMath::logBesselI(k, 2 * sqrtMu1 * sqrtMu2);
    y += 0.5 * k * (logMu1 - logMu2);
    y -= mu1 + mu2;
    return y;
}

template < typename IntType >
double SkellamRand<IntType>::F(const IntType & k) const
{
    return (k < 0) ? RandMath::MarcumP(-k, mu1, mu2, sqrtMu1, sqrtMu2, logMu1, logMu2) : RandMath::MarcumQ(k + 1, mu2, mu1, sqrtMu2, sqrtMu1, logMu2, logMu1);
}

template < typename IntType >
double SkellamRand<IntType>::S(const IntType & k) const
{
    return (k < 0) ? RandMath::MarcumQ(-k, mu1, mu2, sqrtMu1, sqrtMu2, logMu1, logMu2) : RandMath::MarcumP(k + 1, mu2, mu1, sqrtMu2, sqrtMu1, logMu2, logMu1);
}

template < typename IntType >
IntType SkellamRand<IntType>::Variate() const
{
    return X.Variate() - Y.Variate();
}

template < typename IntType >
void SkellamRand<IntType>::Sample(std::vector<IntType> &outputData) const
{
    X.Sample(outputData);
    for (IntType & var : outputData)
        var -= Y.Variate();
}

template < typename IntType >
void SkellamRand<IntType>::Reseed(unsigned long seed) const
{
    X.Reseed(seed);
    Y.Reseed(seed + 1);
}

template < typename IntType >
long double SkellamRand<IntType>::Mean() const
{
    return mu1 - mu2;
}

template < typename IntType >
long double SkellamRand<IntType>::Variance() const
{
    return mu1 + mu2;
}

template < typename IntType >
IntType SkellamRand<IntType>::Median() const
{
    return DiscreteDistribution<IntType>::quantileImpl(0.5, mu1 - mu2);
}

template < typename IntType >
std::complex<double> SkellamRand<IntType>::CFImpl(double t) const
{
    double cosT = std::cos(t), sinT = std::sin(t);
    double x = (cosT - 1) * (mu1 + mu2);
    double y = sinT * (mu1 - mu2);
    std::complex<double> z(x, y);
    return std::exp(z);
}

template < typename IntType >
IntType SkellamRand<IntType>::Mode() const
{
    return Mean();
}

template < typename IntType >
long double SkellamRand<IntType>::Skewness() const
{
    return (mu1 - mu2) / std::pow(mu1 + mu2, 1.5);
}

template < typename IntType >
long double SkellamRand<IntType>::ExcessKurtosis() const
{
    return 1.0 / (mu1 + mu2);
}


template class SkellamRand<int>;
template class SkellamRand<long int>;
template class SkellamRand<long long int>;

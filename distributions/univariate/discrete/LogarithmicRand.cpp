#include "LogarithmicRand.h"
#include "../continuous/UniformRand.h"

template < typename IntType >
LogarithmicRand<IntType>::LogarithmicRand(double probability)
{
    SetProbability(probability);
}

template < typename IntType >
String LogarithmicRand<IntType>::Name() const
{
    return "Logarithmic(" + this->toStringWithPrecision(GetProbability()) + ")";
}

template < typename IntType >
void LogarithmicRand<IntType>::SetProbability(double probability)
{
    if (probability <= 0.0 || probability >= 1.0)
        throw std::invalid_argument("Logarithmic distribution: probability parameter should in interval (0, 1)");
    p = probability;
    logProb = std::log(p);
    log1mProb = std::log1pl(-p);
}

template < typename IntType >
double LogarithmicRand<IntType>::P(const IntType & k) const
{
    return (k < 1) ? 0.0 : -std::pow(p, k) / (k * log1mProb);
}

template < typename IntType >
double LogarithmicRand<IntType>::logP(const IntType & k) const
{
    return (k < 1) ? -INFINITY : k * logProb - std::log(-k * log1mProb);
}

template < typename IntType >
double LogarithmicRand<IntType>::betaFun(IntType a) const
{
    double denom = a + 1;
    double sum = p * p / (a + 2) + p / (a + 1) + 1.0 / a;
    double add = 1;
    int i = 3;
    do {
        add = std::exp(i * logProb) / (++denom);
        sum += add;
        ++i;
    } while (add > MIN_POSITIVE * sum);
    return std::exp(a * logProb) * sum;
}

template < typename IntType >
double LogarithmicRand<IntType>::F(const IntType & k) const
{
    return (k < 1) ? 0.0 : 1 + betaFun(k + 1) / log1mProb;
}

template < typename IntType >
double LogarithmicRand<IntType>::S(const IntType & k) const
{
    return (k < 1) ? 1.0 : -betaFun(k + 1) / log1mProb;
}

template < typename IntType >
IntType LogarithmicRand<IntType>::Variate() const
{
    /// Kemp's second accelerated generator
    /// p. 548, "Non-Uniform Random Variate Generation" by Luc Devroye
    float V = UniformRand<float>::StandardVariate(this->localRandGenerator);
    if (V >= p)
        return 1.0;
    float U = UniformRand<float>::StandardVariate(this->localRandGenerator);
    double y = -std::expm1l(U * log1mProb);
    if (V > y)
        return 1.0;
    if (V > y * y)
        return 2.0;
    return std::floor(1.0 + std::log(V) / std::log(y));
}

template < typename IntType >
long double LogarithmicRand<IntType>::Mean() const
{
    return -p / (1.0 - p) / log1mProb;
}

template < typename IntType >
long double LogarithmicRand<IntType>::Variance() const
{
    long double var = p / log1mProb + 1;
    var /= log1mProb;
    var *= p;
    long double q = 1.0 - p;
    return -var / (q * q);
}

template < typename IntType >
std::complex<double> LogarithmicRand<IntType>::CFImpl(double t) const
{
    std::complex<double> y(std::cos(t), std::sin(t));
    y = std::log(1.0 - p * y);
    return y / log1mProb;
}

template < typename IntType >
IntType LogarithmicRand<IntType>::Mode() const
{
    return 1.0;
}

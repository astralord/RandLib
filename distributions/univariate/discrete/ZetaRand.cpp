#include "ZetaRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ParetoRand.h"

template < typename IntType >
ZetaRand<IntType>::ZetaRand(double exponent)
{
    SetExponent(exponent);
}

template < typename IntType >
String ZetaRand<IntType>::Name() const
{
    return "Zeta(" + this->toStringWithPrecision(GetExponent()) + ")";
}

template < typename IntType >
void ZetaRand<IntType>::SetExponent(double exponent)
{
    if (exponent <= 1.0)
        throw std::invalid_argument("Zeta distribution: exponent should be larger than 1");
    s = exponent;
    sm1 = s - 1.0;
    zetaS = std::riemann_zetal(s);
    logZetaS = std::log(zetaS);
    b = -std::expm1l(-sm1 * M_LN2);
}

template < typename IntType >
double ZetaRand<IntType>::P(const IntType &k) const
{
    return (k < 1) ? 0.0 : std::exp(logP(k));
}

template < typename IntType >
double ZetaRand<IntType>::logP(const IntType & k) const
{
    return (k < 1) ? -INFINITY : -logZetaS - s * std::log(k);
}

template < typename IntType >
double ZetaRand<IntType>::F(const IntType & k) const
{
    return (k < 1) ? 0.0 : RandMath::harmonicNumber(s, k) / zetaS;
}

template < typename IntType >
IntType ZetaRand<IntType>::Variate() const
{
    /// Luc Devroye, p. 551
    /// rejection sampling from rounded down Pareto distribution
    size_t iter = 0;
    do {
        IntType X = std::floor(ParetoRand<float>::StandardVariate(sm1, this->localRandGenerator));
        float T = std::pow(1.0 + 1.0 / X, sm1);
        float V = UniformRand<float>::StandardVariate(this->localRandGenerator);
        /// there was a typo in the book - '<=' instead of '>'
        if (V * X * (T - 1) <= b * T )
            return X;
    } while (++iter <= ProbabilityDistribution<IntType>::MAX_ITER_REJECTION);
    return -1; /// return if algorithm doesn't work
}

template < typename IntType >
long double ZetaRand<IntType>::Mean() const
{
    return (s > 2) ? std::riemann_zetal(sm1) / zetaS : INFINITY;
}

template < typename IntType >
long double ZetaRand<IntType>::Variance() const
{
    if (s <= 3)
        return INFINITY;
    double y = Mean();
    double z = std::riemann_zetal(s - 2) / zetaS;
    return z - y * y;
}

template < typename IntType >
IntType ZetaRand<IntType>::Mode() const
{
    return 1;
}

template < typename IntType >
long double ZetaRand<IntType>::Skewness() const
{
    if (s <= 4)
        return INFINITY;
    long double z1 = std::riemann_zetal(sm1), z1Sq = z1 * z1;
    long double z2 = std::riemann_zetal(s - 2);
    long double z3 = std::riemann_zetal(s - 3);
    long double z = zetaS, zSq = z * z;
    long double numerator = zSq * z3;
    numerator -= 3 * z2 * z1 * z;
    numerator += 2 * z1 * z1Sq;
    long double denominator = z * z2 - z1Sq;
    denominator = std::pow(denominator, 1.5);
    denominator *= zSq;
    return numerator / denominator;
}

template < typename IntType >
long double ZetaRand<IntType>::ExcessKurtosis() const
{
    if (s <= 5)
        return INFINITY;
    long double mean = Mean();
    long double secondMoment = this->SecondMoment();
    long double thirdMoment = ThirdMoment();
    long double fourthMoment = FourthMoment();
    long double meanSq = mean * mean;
    long double variance = secondMoment - meanSq;
    long double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    long double denominator = variance * variance;
    return numerator / denominator - 3.0;
}

template < typename IntType >
long double ZetaRand<IntType>::Moment(int n) const
{
    return (s > n + 1) ? std::riemann_zetal(s - n) / zetaS : INFINITY;
}


#include "ZetaRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ParetoRand.h"

namespace
{

long double zetaTail(double exponent, long double firstTerm)
{
    static constexpr long double B2jOverFact[] = {
        1.0L / 12, -1.0L / 720, 1.0L / 30240, -1.0L / 1209600, 1.0L / 47900160,
        -691.0L / 1307674368000.0L, 1.0L / 74724249600.0L, -3617.0L / 10670622842880000.0L,
        43867.0L / 5109094217170944000.0L, -174611.0L / 802857662698291200000.0L
    };
    long double nPow = std::pow(firstTerm, -static_cast<long double>(exponent));
    long double sum = firstTerm * nPow / (exponent - 1.0L) + 0.5L * nPow;
    long double rising = exponent, term = nPow / firstTerm;
    for (int j = 0; j < 10; ++j) {
        sum += B2jOverFact[j] * rising * term;
        rising *= (exponent + 2 * j + 1) * (exponent + 2 * j + 2);
        term /= firstTerm * firstTerm;
    }
    return sum;
}

}

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
    if (!(exponent > 1.0))
        throw std::invalid_argument("Zeta distribution: exponent should be greater than 1, but it's equal to " + std::to_string(exponent));
    s = exponent;
    sm1 = s - 1.0;
    zetaS = RandMath::riemannZeta(s);
    logZetaS = std::log(zetaS);
    b = -std::expm1l(-sm1 * randlib::consts::ln2);
}

template < typename IntType >
double ZetaRand<IntType>::logP(const IntType & k) const
{
    return (k < 1) ? -INFINITY : -logZetaS - s * std::log(k);
}

template < typename IntType >
double ZetaRand<IntType>::F(const IntType & k) const
{
    if (k < 1)
        return 0.0;
    return (k <= 100000) ? RandMath::harmonicNumber(s, k) / zetaS : 1.0 - S(k);
}

template < typename IntType >
double ZetaRand<IntType>::S(const IntType & k) const
{
    if (k < 1)
        return 1.0;
    if (k <= 100000)
        return 1.0 - F(k);
    return static_cast<double>(zetaTail(s, static_cast<long double>(k) + 1.0L) / zetaS);
}

template < typename IntType >
IntType ZetaRand<IntType>::Variate() const
{
    /// Luc Devroye, p. 551
    /// rejection sampling from rounded down Pareto distribution
    size_t iter = 0;
    do {
        float X = std::floor(ParetoRand<float>::StandardVariate(sm1, this->localRandGenerator));
        float V = UniformRand<float>::StandardVariate(this->localRandGenerator);
        if (X < 1e4) {
            float T = std::pow(1.0 + 1.0 / X, sm1);
            /// there was a typo in the book - '<=' instead of '>'
            if (V * X * (T - 1) <= b * T)
                return X;
        }
        else {
            long double logT = sm1 * std::log1pl(1.0 / X);
            long double TM1 = std::expm1l(logT);
            if (V * X * TM1 <= b * std::exp(logT))
                return X < this->MaxValue() ? (IntType)X : this->MaxValue();
        }

    } while (++iter <= ProbabilityDistribution<IntType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Zeta distribution: sampling failed");
}

template < typename IntType >
long double ZetaRand<IntType>::Mean() const
{
    return (s > 2) ? RandMath::riemannZeta(sm1) / zetaS : INFINITY;
}

template < typename IntType >
long double ZetaRand<IntType>::Variance() const
{
    if (s <= 3)
        return INFINITY;
    double y = Mean();
    double z = RandMath::riemannZeta(s - 2) / zetaS;
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
    long double z1 = RandMath::riemannZeta(sm1), z1Sq = z1 * z1;
    long double z2 = RandMath::riemannZeta(s - 2);
    long double z3 = RandMath::riemannZeta(s - 3);
    long double z = zetaS, zSq = z * z;
    long double logskew = zSq * z3;
    logskew -= 3 * z2 * z1 * z;
    logskew += 2 * z1 * z1Sq;
    logskew = std::log(logskew);
    logskew -= 1.5 * std::log(z * z2 - z1Sq);
    logskew -= 2 * logZetaS;
    return std::exp(logskew);
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
    return (s > n + 1) ? RandMath::riemannZeta(s - n) / zetaS : INFINITY;
}

template class ZetaRand<int>;
template class ZetaRand<long int>;
template class ZetaRand<long long int>;

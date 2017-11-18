#include "ZetaRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ParetoRand.h"

ZetaRand::ZetaRand(double exponent)
{
    SetExponent(exponent);
}

String ZetaRand::Name() const
{
    return "Zeta(" + toStringWithPrecision(GetExponent()) + ")";
}

void ZetaRand::SetExponent(double exponent)
{
    if (exponent <= 1.0)
        throw std::invalid_argument("Zeta distribution: exponent should be larger than 1");
    s = exponent;
    sm1 = s - 1.0;
    zetaS = std::riemann_zetal(s);
    logZetaS = std::log(zetaS);
    b = -std::expm1l(-sm1 * M_LN2);
}

double ZetaRand::P(const int & k) const
{
    return (k < 1) ? 0.0 : std::exp(logP(k));
}

double ZetaRand::logP(const int & k) const
{
    return (k < 1) ? -INFINITY : -logZetaS - s * std::log(k);
}

double ZetaRand::F(const int & k) const
{
    return (k < 1) ? 0.0 : RandMath::harmonicNumber(s, k) / zetaS;
}

int ZetaRand::Variate() const
{
    /// Luc Devroye, p. 551
    /// rejection sampling from rounded down Pareto distribution
    int iter = 0;
    do {
        double X = std::floor(ParetoRand::StandardVariate(sm1, localRandGenerator));
        double T = std::pow(1.0 + 1.0 / X, sm1);
        double V = UniformRand::StandardVariate(localRandGenerator);
        /// there was a typo in the book - '<=' instead of '>'
        if (V * X * (T - 1) <= b * T )
            return X;
    } while (++iter <= MAX_ITER_REJECTION);
    return -1; /// return if algorithm doesn't work
}

long double ZetaRand::Mean() const
{
    return (s > 2) ? std::riemann_zetal(sm1) / zetaS : INFINITY;
}

long double ZetaRand::Variance() const
{
    if (s <= 3)
        return INFINITY;
    double y = Mean();
    double z = std::riemann_zetal(s - 2) / zetaS;
    return z - y * y;
}

int ZetaRand::Mode() const
{
    return 1;
}

long double ZetaRand::Skewness() const
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

long double ZetaRand::ExcessKurtosis() const
{
    if (s <= 5)
        return INFINITY;
    long double mean = Mean();
    long double secondMoment = SecondMoment();
    long double thirdMoment = ThirdMoment();
    long double fourthMoment = FourthMoment();
    long double meanSq = mean * mean;
    long double variance = secondMoment - meanSq;
    long double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    long double denominator = variance * variance;
    return numerator / denominator - 3.0;
}

long double ZetaRand::Moment(int n) const
{
    return (s > n + 1) ? std::riemann_zetal(s - n) / zetaS : INFINITY;
}


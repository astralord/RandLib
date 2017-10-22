#include "ZetaRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ParetoRand.h"

ZetaRand::ZetaRand(double exponent)
{
    SetExponent(exponent);
}

std::string ZetaRand::Name() const
{
    return "Zeta(" + toStringWithPrecision(GetExponent()) + ")";
}

void ZetaRand::SetExponent(double exponent)
{
    s = exponent <= 1.0 ? 2.0 : exponent;
    sm1 = s - 1.0;
    zetaSInv = 1.0 / RandMath::zetaRiemann(s);
    b = -std::expm1(-sm1 * M_LN2);
}

double ZetaRand::P(const int & k) const
{
    return (k < 1) ? 0.0 : zetaSInv / std::pow(k, s);
}

double ZetaRand::logP(const int & k) const
{
    return (k < 1) ? -INFINITY : std::log(zetaSInv) - s * std::log(k); // log(1/zeta(s)) can be hashed
}

double ZetaRand::F(const int & k) const
{
    return (k < 1) ? 0.0 : zetaSInv * RandMath::harmonicNumber(s, k);
}

int ZetaRand::Variate() const
{
    /// Luc Devroye, p. 551
    /// rejection sampling from rounded down Pareto distribution
    int iter = 0;
    do {
        double X = std::floor(ParetoRand::StandardVariate(sm1));
        double T = std::pow(1.0 + 1.0 / X, sm1);
        double V = UniformRand::StandardVariate();
        /// there was a typo in the book - '<=' instead of '>'
        if (V * X * (T - 1) <= b * T )
            return X;
    } while (++iter <= MAX_ITER_REJECTION);
    return -1; /// return if algorithm doesn't work
}

double ZetaRand::Mean() const
{
    return (s > 2) ? zetaSInv * RandMath::zetaRiemann(sm1) : INFINITY;
}

double ZetaRand::Variance() const
{
    if (s <= 3)
        return INFINITY;
    double y = Mean();
    double z = zetaSInv * RandMath::zetaRiemann(s - 2);
    return z - y * y;
}

int ZetaRand::Mode() const
{
    return 1;
}

double ZetaRand::Skewness() const
{
    if (s <= 4)
        return INFINITY;
    double z1 = RandMath::zetaRiemann(sm1), z1Sq = z1 * z1;
    double z2 = RandMath::zetaRiemann(s - 2);
    double z3 = RandMath::zetaRiemann(s - 3);
    double z = 1.0 / zetaSInv, zSq = z * z;
    double numerator = zSq * z3;
    numerator -= 3 * z2 * z1 * z;
    numerator += 2 * z1 * z1Sq;
    double denominator = z * z2 - z1Sq;
    denominator = std::pow(denominator, 1.5);
    denominator *= zSq;
    return numerator / denominator;
}

double ZetaRand::ExcessKurtosis() const
{
    if (s <= 5)
        return INFINITY;
    double mean = Mean();
    double secondMoment = SecondMoment();
    double thirdMoment = ThirdMoment();
    double fourthMoment = FourthMoment();
    double meanSq = mean * mean;
    double variance = secondMoment - meanSq;
    double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    double denominator = variance * variance;
    return numerator / denominator - 3.0;
}

double ZetaRand::Moment(int n) const
{
    return (s > n + 1) ? zetaSInv * RandMath::zetaRiemann(s - n) : INFINITY;
}


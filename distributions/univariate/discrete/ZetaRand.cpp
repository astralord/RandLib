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
    b = 1.0 - std::pow(2.0, -sm1);
}

double ZetaRand::P(int k) const
{
    if (k < 1)
        return 0.0;
    return zetaSInv / std::pow(k, s);
}

double ZetaRand::F(int k) const
{
    if (k < 1.0)
        return 0.0;
    return zetaSInv * RandMath::harmonicNumber(s, k);
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
    if (s <= 2)
        return INFINITY;
    return zetaSInv * RandMath::zetaRiemann(sm1);
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
    return DiscreteDistribution::ExcessKurtosis(); // TODO: implement
}


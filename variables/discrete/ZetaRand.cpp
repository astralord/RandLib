#include "ZetaRand.h"
#include "../continuous/UniformRand.h"

ZetaRand::ZetaRand(double exponent)
{
    setExponent(exponent);
}

std::string ZetaRand::name()
{
    return "Zeta(" + toStringWithPrecision(getExponent()) + ")";
}

void ZetaRand::setExponent(double exponent)
{
    s = exponent;
    /// sanity check
    if (s <= 1.0)
        s = 2.0; /// default value
    zetaSInv = 1.0 / RandMath::zetaRiemann(s);
    b = std::pow(2.0, s - 1.0);
}

double ZetaRand::P(int k) const
{
    if (k < 1)
        return 0.0;
    return zetaSInv / std::pow(k, s);
}

double ZetaRand::F(double x) const
{
    if (x < 1.0)
        return 0.0;
    int k = static_cast<int>(std::floor(x));
    return zetaSInv * RandMath::harmonicNumber(s, k);
}

double ZetaRand::variate() const
{
    /// Luc Devroye, p. 551
    int iter = 0;
    double sm1 = s - 1.0;
    do {
        /// faster than calling Pareto
        double U = UniformRand::standardVariate();
        double X = std::pow(U, -1.0 / sm1);
        X = std::floor(X);
        double T = std::pow(1.0 + 1.0 / X, sm1);
        double V = UniformRand::standardVariate();
        /// there was typo in the book - '<=' instead of '>'
        if (V * X * (T - 1) * b <= T * (b - 1))
            return X;
    } while (++iter < 1e9);
    return NAN; /// doesn't work
}

double ZetaRand::Mean() const
{
    if (s <= 2)
        return INFINITY;
    return zetaSInv * RandMath::zetaRiemann(s - 1);
}

double ZetaRand::Variance() const
{
    if (s <= 3)
        return INFINITY;
    double y = Mean();
    double z = zetaSInv * RandMath::zetaRiemann(s - 2);
    return z - y * y;
}

double ZetaRand::Mode() const
{
    return 1.0;
}


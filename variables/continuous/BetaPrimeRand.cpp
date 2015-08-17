#include "BetaPrimeRand.h"

BetaPrimeRand::BetaPrimeRand(double shape1, double shape2)
    : BetaRand(shape1, shape2)
{
    setName();
}

void BetaPrimeRand::setName()
{
    nameStr = "Beta Prime(" + toStringWithPrecision(getAlpha()) + ", " + toStringWithPrecision(getBeta()) + ")";
}

double BetaPrimeRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double rv = std::pow(x, alpha - 1);
    rv *= std::pow(1 + x, -alpha - beta);
    return pdfCoef * rv;
}

double BetaPrimeRand::F(double x) const
{
    //TODO
    return x;
}

double BetaPrimeRand::variate() const
{
    double b = BetaRand::variate();
    return b / (1 - b);
}

#include "WeibullRand.h"

WeibullRand::WeibullRand(double scale, double shape)
{
    setParameters(scale, shape);
}

void WeibullRand::setName()
{
    nameStr = "Weibull(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void WeibullRand::setParameters(double scale, double shape)
{
    l = std::max(scale, MIN_POSITIVE);
    k = std::max(shape, MIN_POSITIVE);
    lInv = 1.0 / l;

    setName();
}

double WeibullRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double xAdj = x * lInv;
    double xAdjPow = std::pow(xAdj, k - 1);
    return k * lInv * xAdjPow * std::exp(-xAdj * xAdjPow);
}

double WeibullRand::F(double x) const
{
    if (x <= 0)
        return 0;
    return 1 - std::exp(-std::pow(x * lInv, k));
}

double WeibullRand::variate() const
{
    return l * std::pow(ExponentialRand::standardVariate(), 1.0 / k);
}

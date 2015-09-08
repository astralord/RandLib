#include "YuleRand.h"

YuleRand::YuleRand(double shape)
{
    setShape(shape);
}

std::string YuleRand::name()
{
    return "Yule(" + toStringWithPrecision(getShape()) + ")";
}

void YuleRand::setShape(double shape)
{
    ro = std::max(shape, MIN_POSITIVE);
    gamma1pRo = std::tgamma(ro + 1);
}

double YuleRand::P(int k) const
{
    if (k < 1)
        return 0;
    return ro * gamma1pRo * std::tgamma(k) / std::tgamma(k + ro + 1);
}

double YuleRand::F(double x) const
{
    if (x < 1)
        return 0;
    double k = std::floor(x);
    return 1 - k * gamma1pRo * std::tgamma(k) / std::tgamma(k + ro + 1);
}

double YuleRand::variate() const
{
    return YuleRand::variate(ro);
}

double YuleRand::variate(double shape)
{
    double prob = 1.0 / ParetoRand::variate(shape, 1.0);
    return GeometricRand::variate(prob) + 1;
}


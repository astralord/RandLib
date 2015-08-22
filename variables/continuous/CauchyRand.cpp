#include "CauchyRand.h"

CauchyRand::CauchyRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

void CauchyRand::setName()
{
    nameStr = "Cauchy(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void CauchyRand::setLocation(double location)
{
    x0 = location;
    setName();
}

void CauchyRand::setScale(double scale)
{
    gamma = std::max(scale, MIN_POSITIVE);
    gammaInv = 1.0 / gamma;
    setName();
}

double CauchyRand::f(double x) const
{
    double y = x - x0;
    y *= y;
    y *= gammaInv;
    y += gamma;
    return M_1_PI / y;
}

double CauchyRand::F(double x) const
{
    double y = x - x0;
    y *= gammaInv;
    y = std::atan(y);
    y *= M_1_PI;
    return y + .5;
}

double CauchyRand::variate() const
{
    return x0 + gamma * standardVariate();
}

double CauchyRand::variate(double location, double scale)
{
    return location + scale * standardVariate();
}

double CauchyRand::standardVariate()
{
    double x, y;
    do {
        x = UniformRand::variate(-1, 1);
        y = UniformRand::variate(-1, 1);
    } while (x * x + y * y > 1.0 || y == 0.0);
    return x / y;
}

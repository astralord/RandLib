#include "CauchyRand.h"

CauchyRand::CauchyRand(double location, double scale) :
    U(-1, 1)
{
    setLocation(location);
    setScale(scale);
}

void CauchyRand::setLocation(double location)
{
    x0 = location;
}

void CauchyRand::setScale(double scale)
{
    gamma = std::max(scale, MIN_POSITIVE);
    gammaInv = 1.0 / gamma;
}

double CauchyRand::f(double x) const
{
    double y = x - x0; /// x - x0
    y *= y; /// (x - x0) ^ 2
    y *= gammaInv; /// (x - x0) ^ 2 / gamma
    y += gamma; /// gamma + (x - x0) ^ 2 / gamma
    return M_1_PI / y; /// gamma / (pi * (gamma ^ 2 + (x - x0) ^ 2))
}

double CauchyRand::F(double x) const
{
    double y = x - x0; /// x - x0
    y *= gammaInv; /// (x - x0) / gamma
    y = std::atan(y); /// atan((x - x0) / gamma)
    y *= M_1_PI; /// atan((x - x0) / gamma) / pi
    return y + .5; /// atan((x - x0) / gamma) / pi + 0.5
}

double CauchyRand::value()
{
    double u1, u2;
    do {
        u1 = U.value();
        u2 = U.value();
    } while (u1 * u1 + u2 * u2 > 1.0 || u2 == 0.0);
    u1 /= u2;
    u1 *= gamma;
    u1 += x0;
    return u1;
}

#include "WaldRand.h"

WaldRand::WaldRand(double mean, double shape)
{
    setParameters(mean, shape);
}

void WaldRand::setParameters(double mean, double shape)
{
    mu = std::max(mean, MIN_POSITIVE);
    l = std::max(shape, MIN_POSITIVE);
    pdfCoef = std::sqrt(.5 * l * M_1_PI);
    cdfCoef = std::exp(2 * l / mu);
}

double WaldRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double xInv = 1.0 / x;
    double y = xInv * std::sqrt(xInv);
    double z = (x - mu);
    z *= z;
    z *= -.5 * l * xInv / (mu * mu);
    z = std::exp(z);
    return pdfCoef * y * z;
}

double WaldRand::F(double x) const
{
    if (x <= 0)
        return 0;
    double a = std::sqrt(l / x);
    double b = x / mu;
    double y = X.F(a * b - a);
    double z = X.F(-a * b - a);
    return y + cdfCoef * z;
}

double WaldRand::value()
{
    double y = X.value();
    y *= y;
    double my = mu * y;
    double x = 4 * l + my;
    x *= my;
    x = std::sqrt(x);
    x = my - x;
    x *= .5 / l;
    ++x;
    if (U.value() <= 1.0 / (1 + x))
        return mu * x;
    return mu / x;
}

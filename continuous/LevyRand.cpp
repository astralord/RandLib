#include "LevyRand.h"

LevyRand::LevyRand(double location, double scale) :
    X(0.0, 1.0)
{
    setLocation(location);
    setScale(scale);
}

void LevyRand::setLocation(double location)
{
    mu = location;
}

void LevyRand::setScale(double scale)
{
    c_2 = scale > 0.0 ? scale : 0.0;
    sqrtc_2pi = std::sqrt(c_2);
    X.setSigma(1.0 / sqrtc_2pi); /// X ~ N(0, 1 / c ^ (0.5))
    sqrtc_2pi *= M_1_SQRT2PI;
    c_2 *= .5;
}

double LevyRand::pdf(double x) const
{
    if (x <= mu)
        return 0;
    double xInv = 1.0 / (x - mu);
    double y = -c_2 * xInv;
    y = std::exp(y);
    y *= xInv;
    y *= std::sqrt(xInv);
    return sqrtc_2pi * y;
}

double LevyRand::cdf(double x) const
{
    if (x <= mu)
        return 0;
    double y = x - mu;
    y = c_2 / y;
    y = std::sqrt(y);
    return erfc(y);
}

double LevyRand::value()
{
    double rv = X.value();
    rv *= rv;
    rv = 1.0 / rv;
    return mu + rv;
}


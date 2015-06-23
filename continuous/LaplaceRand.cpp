#include "LaplaceRand.h"

LaplaceRand::LaplaceRand(double location, double scale) :
    X(1)
{
    setLocation(location);
    setScale(scale);
}

void LaplaceRand::setLocation(double location)
{
    mu = location;
}

void LaplaceRand::setScale(double scale)
{
    b = std::max(scale, MIN_POSITIVE);
    bInv = 1.0 / b;
    X.setRate(bInv);
}

double LaplaceRand::f(double x) const
{
    double y = -std::fabs(x - mu);
    y *= bInv;
    y = std::exp(y);
    y *= bInv;
    return .5 * y;
}

double LaplaceRand::F(double x) const
{
    double y = x - mu;
    y *= bInv;
    if (x < mu)
        return .5 * std::exp(y);
    y = -.5 * std::exp(-y);
    return y + 1;
}

double LaplaceRand::value()
{
    double rv = X.value(); /// ~ Exp(1/b)
    rv -= X.value(); /// ~ Laplace(0, b)
    return mu + rv;
}

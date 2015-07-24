#include "LaplaceRand.h"

LaplaceRand::LaplaceRand(double location, double scale)
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
    double e = X.value();
    return mu + (((signed)B.getRand() > 0) ? e : -e);
}

#include "LogisticRand.h"

LogisticRand::LogisticRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

void LogisticRand::setLocation(double location)
{
    mu = location;
}

void LogisticRand::setScale(double scale)
{
    s = scale;
}

double LogisticRand::f(double x) const
{
    double numen = std::exp((mu - x) / s);
    double denom = (1 + numen);
    denom *= denom;
    denom *= s;
    return numen / denom;
}

double LogisticRand::F(double x) const
{
    double expX = std::exp((mu - x) / s);
    return 1.0 / (1 + expX);
}

double LogisticRand::value()
{
    double rv = U.value();
    rv = rv / (1 - rv);
    rv = std::log(rv);
    return mu + s * rv;
}

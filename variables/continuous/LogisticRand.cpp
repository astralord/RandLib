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
    double numerator = std::exp((mu - x) / s);
    double denominator = (1 + numerator);
    denominator *= denominator;
    denominator *= s;
    return numerator / denominator;
}

double LogisticRand::F(double x) const
{
    double expX = std::exp((mu - x) / s);
    return 1.0 / (1 + expX);
}

double LogisticRand::variate() const
{
    return mu + s * std::log(1.0 / UniformRand::standardVariate() - 1);
}

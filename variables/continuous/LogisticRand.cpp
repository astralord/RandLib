#include "LogisticRand.h"

LogisticRand::LogisticRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string LogisticRand::name()
{
    return "Logistic(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LogisticRand::setLocation(double location)
{
    mu = location;
}

void LogisticRand::setScale(double scale)
{
    s = scale;
    if (s <= 0)
        s = MIN_POSITIVE;
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

double LogisticRand::Variance() const
{
    double sPi = s * M_PI;
    return sPi * sPi / 3;
}

double LogisticRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return -INFINITY;
    if (p == 1)
        return INFINITY;
    return mu - s * std::log(1.0 / p - 1);
}

double LogisticRand::Median() const
{
    return mu;
}

double LogisticRand::Mode() const
{
    return mu;
}

double LogisticRand::Skewness() const
{
    return 0;
}

double LogisticRand::ExcessKurtosis() const
{
    return 1.2;
}

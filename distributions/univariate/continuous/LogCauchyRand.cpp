#include "LogCauchyRand.h"

LogCauchyRand::LogCauchyRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string LogCauchyRand::name()
{
    return "Log-Cauchy(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LogCauchyRand::setLocation(double location)
{
    X.setLocation(location);
    expMu = std::exp(X.Mean());
}

void LogCauchyRand::setScale(double scale)
{
    X.setScale(scale);
    expVar = std::exp(X.Variance());
}

double LogCauchyRand::f(double x) const
{
    return (x > 0) ? X.f(std::log(x)) / x : 0.0;
}

double LogCauchyRand::F(double x) const
{
    return (x > 0) ? X.F(std::log(x)) : 0.0;
}

double LogCauchyRand::variate() const
{
    return std::exp(X.variate());
}

double LogCauchyRand::Mean() const
{
    return NAN;
}

double LogCauchyRand::Variance() const
{
    return INFINITY;
}

double LogCauchyRand::Quantile(double p) const
{
    return (p == 0) ? 0.0 : std::exp(X.Quantile(p));
}

double LogCauchyRand::Median() const
{
    return expMu;
}

double LogCauchyRand::Mode() const
{
    return 0.0;
}

double LogCauchyRand::Skewness() const
{
    return NAN;
}

double LogCauchyRand::ExcessKurtosis() const
{
    return NAN;
}

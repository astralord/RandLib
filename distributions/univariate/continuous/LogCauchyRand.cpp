#include "LogCauchyRand.h"

LogCauchyRand::LogCauchyRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string LogCauchyRand::name() const
{
    return "Log-Cauchy(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LogCauchyRand::setLocation(double location)
{
    X.setLocation(location);
    expMu = std::exp(X.getLocation());
}

void LogCauchyRand::setScale(double scale)
{
    X.setScale(scale);
}

double LogCauchyRand::f(double x) const
{
    if (x == 0)
        return INFINITY;
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

double LogCauchyRand::quantileImpl(double p) const
{
    return std::exp(X.Quantile(p));
}

double LogCauchyRand::quantileImpl1m(double p) const
{
    return std::exp(X.Quantile1m(p));
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

#include "StableProcess.h"

StableProcess::StableProcess(double exponent, double skewness, double scale, double location, double deltaT) :
    StochasticProcess(deltaT),
    X(exponent, skewness, scale, location),
    dtCoef(std::pow(dt, X.getInvExponent()))
{

}

void StableProcess::nextImpl()
{
    currentValue += dtCoef * X.variate();
}

void StableProcess::nextImpl(double deltaT)
{
    currentValue += std::pow(deltaT, X.getInvExponent()) * X.variate();
}

double StableProcess::MeanImpl(double t) const
{
    double alpha = X.getExponent();
    if (alpha > 1)
        return currentValue + std::pow(t - currentTime, X.getInvExponent()) * X.Mean();
    return (alpha == 1) ? NAN : INFINITY;
}

double StableProcess::VarianceImpl(double t) const
{
    return (X.getExponent() == 2.0) ? t - currentTime : INFINITY;
}

double StableProcess::QuantileImpl(double t, double p) const
{
    return currentValue + std::pow(t - currentTime, X.getInvExponent()) * X.Quantile(p);
}

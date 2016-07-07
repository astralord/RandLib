#include "StableProcess.h"

StableProcess::StableProcess(double exponent, double skewness, double scale, double location, double deltaT) :
    StochasticProcess(deltaT),
    X(exponent, skewness, 1.0, 0.0),
    dtCoef(std::pow(dt, X.getInvExponent())),
    drift(location),
    volatility(scale)
{

}

void StableProcess::nextImpl()
{
    currentValue += drift * dt + dtCoef * volatility * X.variate();
}

void StableProcess::nextImpl(double deltaT)
{
    currentValue += drift * dt + volatility * std::pow(deltaT, X.getInvExponent()) * X.variate();
}

double StableProcess::MeanImpl(double t) const
{
    return (X.getExponent() > 1) ? currentValue + drift * (t - currentTime) : NAN;
}

double StableProcess::VarianceImpl(double t) const
{
    return (X.getExponent() == 2.0) ? volatility * volatility * (t - currentTime) : INFINITY;
}

double StableProcess::QuantileImpl(double t, double p) const
{
    double q = X.Quantile(p);
    q *= volatility;
    q *= std::pow(t - currentTime, X.getInvExponent());
    q += drift * (t - currentTime);
    return currentValue + q;
}

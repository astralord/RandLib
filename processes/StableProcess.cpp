#include "StableProcess.h"

StableProcess::StableProcess(double exponent, double skewness, double drift, double volatility, double deltaT) :
    StochasticProcess(deltaT),
    X(exponent, skewness, 1.0, 0.0),
    dtCoef(std::pow(dt, X.getInvExponent())),
    mu(drift),
    sigma(volatility)
{

}

void StableProcess::nextImpl()
{
    currentValue += mu * dt + dtCoef * sigma * X.variate();
}

void StableProcess::nextImpl(double deltaT)
{
    currentValue += mu * dt + sigma * std::pow(deltaT, X.getInvExponent()) * X.variate();
}

double StableProcess::MeanImpl(double t) const
{
    return (X.getExponent() > 1) ? currentValue + mu * (t - currentTime) : NAN;
}

double StableProcess::VarianceImpl(double t) const
{
    return (X.getExponent() == 2.0) ? sigma * sigma * (t - currentTime) : INFINITY;
}

double StableProcess::QuantileImpl(double t, double p) const
{
    double q = X.Quantile(p);
    q *= sigma;
    q *= std::pow(t - currentTime, X.getInvExponent());
    q += mu * (t - currentTime);
    return currentValue + q;
}

void StableProcess::QuantileImpl(const std::vector<double> &t, std::vector<double> &outputData, double p) const
{
    double q = sigma * X.Quantile(p);
    double alphaInv = X.getInvExponent();
    for (size_t i = 0; i != t.size(); ++i)
    {
        outputData[i] = q;
        outputData[i] *= std::pow(t[i] - currentTime, alphaInv);
        outputData[i] += mu * (t[i] - currentTime);
        outputData[i] += currentValue;
    }
}

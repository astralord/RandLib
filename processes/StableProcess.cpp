#include "StableProcess.h"

StableProcess::StableProcess(double exponent, double skewness, double drift, double volatility, double deltaT) :
    StochasticProcess(deltaT),
    X(exponent, skewness, 1.0, 0.0),
    dtCoef(std::pow(dt, 1.0 / X.GetExponent())),
    mu(drift),
    sigma(volatility)
{

}

void StableProcess::NextImpl()
{
    currentValue += mu * dt + dtCoef * sigma * X.Variate();
}

double StableProcess::MeanImpl(double t) const
{
    return (X.GetExponent() > 1) ? currentValue + mu * (t - currentTime) : NAN;
}

double StableProcess::VarianceImpl(double t) const
{
    return (X.GetExponent() == 2.0) ? sigma * sigma * (t - currentTime) : INFINITY;
}

double StableProcess::Quantile(double t, double p) const
{
    if (p < 0 || p > 1 || t < currentTime)
        return NAN;
    double q = X.Quantile(p);
    q *= sigma;
    q *= std::pow(t - currentTime, 1.0 / X.GetExponent());
    q += mu * (t - currentTime);
    return currentValue + q;
}

void StableProcess::Quantile(const std::vector<double> &t, std::vector<double> &outputData, double p) const
{
    /// Assumed that t is sorted
    if (p < 0 || p > 1 || t[0] < currentTime)
        return;
    if (t.size() > outputData.size())
        return;

    double q = sigma * X.Quantile(p);
    double alpha = X.GetExponent();
    for (size_t i = 0; i != t.size(); ++i)
    {
        outputData[i] = q;
        outputData[i] *= std::pow(t[i] - currentTime, 1.0 / alpha);
        outputData[i] += mu * (t[i] - currentTime);
        outputData[i] += currentValue;
    }
}

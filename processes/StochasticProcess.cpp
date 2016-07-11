#include "StochasticProcess.h"

StochasticProcess::StochasticProcess(double deltaT, double initialValue) :
    currentTime(0.0),
    currentValue(initialValue),
    dt(std::max(deltaT, 0.0))
{
}

void StochasticProcess::reset(double initialValue)
{
    currentTime = 0.0;
    currentValue = initialValue;
}

double StochasticProcess::next()
{
    currentTime += dt;
    nextImpl();
    return currentValue;
}

double StochasticProcess::Mean(double t) const
{
    if (t < currentTime)
        return NAN; /// we consider only future time
    return (t > currentTime) ? MeanImpl(t) : currentValue;
}

double StochasticProcess::Variance(double t) const
{
    if (t < currentTime)
        return NAN; /// we consider only future time
    return (t > currentTime) ? VarianceImpl(t) : 0.0;
}

double StochasticProcess::Quantile(double t, double p) const
{
    if (p < 0 || p > 1 || t < currentTime)
        return NAN;
    return QuantileImpl(t, p);
}

void StochasticProcess::Quantile(const std::vector<double> &t, std::vector<double> &outputData, double p) const
{
    /// Assumed that t is sorted
    if (p < 0 || p > 1 || t[0] < currentTime)
        return;
    if (t.size() > outputData.size())
        return;
    return this->QuantileImpl(t, outputData, p);
}

void StochasticProcess::QuantileImpl(const std::vector<double> &t, std::vector<double> &outputData, double p) const
{
    for (size_t i = 0; i != t.size(); ++i)
        outputData[i] = QuantileImpl(t[i], p);
}

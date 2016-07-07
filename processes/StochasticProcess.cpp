#include "StochasticProcess.h"

StochasticProcess::StochasticProcess(double deltaT) :
    currentTime(0.0),
    currentValue(0.0),
    dt(std::max(deltaT, 0.0))
{
}

double StochasticProcess::next()
{
    currentTime += dt;
    nextImpl();
    return currentValue;
}

double StochasticProcess::next(double deltaT)
{
    if (deltaT < 0)
        return NAN;
    if (deltaT > 0) {
        currentTime += deltaT;
        nextImpl(deltaT);
    }
    return currentValue;
}

double StochasticProcess::Mean(double t) const
{
    if (t < currentTime)
        return NAN; /// we consider only future time
    if (t == currentTime)
        return 0.0;
    return MeanImpl(t);
}

double StochasticProcess::Variance(double t) const
{
    if (t < currentTime)
        return NAN; /// we consider only future time
    if (t == currentTime)
        return 0.0;
    return VarianceImpl(t);
}

double StochasticProcess::Quantile(double t, double p) const
{
    if (p < 0 || p > 1 || t <= currentTime)
        return NAN;
    return QuantileImpl(t, p);
}

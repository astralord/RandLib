#include "HestonProcess.h"


HestonProcess::HestonProcess(double drift, double volatilityDrift, double reversionSpeed, double volatility, double initialValue, double volatilityInitialValue, double correlation, double deltaT) :
    StochasticProcess(deltaT, initialValue),
    mu(drift),
    ro(correlation),
    V(volatilityDrift, reversionSpeed, volatility, volatilityInitialValue)
{

}

void HestonProcess::nextImpl()
{
    double v = V.next();
    double rootVT = std::sqrt(v * dt);
    currentValue *= LogNormalRand::variate((mu - 0.5 * v) * dt, rootVT);
}

double HestonProcess::MeanImpl(double t) const
{
    // TODO:
    return t;
}

double HestonProcess::VarianceImpl(double t) const
{
    // TODO:
    return t;
}

double HestonProcess::QuantileImpl(double t, double p) const
{
    // not todo - remove quantile functions
    return t + p;
}

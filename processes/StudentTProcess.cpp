#include "StudentTProcess.h"

StudentTProcess::StudentTProcess(int degree, double drift, double volatility, double deltaT) :
    StochasticProcess(deltaT),
    X(degree)
{
    mu = drift;
    sigma = volatility > 0 ? volatility : 1.0;
}

void StudentTProcess::NextImpl()
{
    // TODO: add dtCoef, so for v = inf it would be Brownian and for v == 1 it would be Cauchy
    currentValue += mu * dt + sigma * X.Variate();
}

double StudentTProcess::MeanImpl(double t) const
{
    return currentValue + mu * (t - currentTime);
}

double StudentTProcess::VarianceImpl(double t) const
{
    return sigma * sigma * (t - currentTime);
}

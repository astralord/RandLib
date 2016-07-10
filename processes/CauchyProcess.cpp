#include "CauchyProcess.h"
#include "../distributions/univariate/continuous/CauchyRand.h"

CauchyProcess::CauchyProcess(double drift, double volatility, double deltaT) :
    StableProcess(1, 0, drift, volatility, deltaT)
{
}

void CauchyProcess::nextImpl()
{
    currentValue += CauchyRand::variate(mu, sigma) * dt;
}

void CauchyProcess::nextImpl(double deltaT)
{
    currentValue += CauchyRand::variate(mu, sigma) * deltaT;
}

double CauchyProcess::QuantileImpl(double t, double p) const
{
    return CauchyRand::quantile(p, currentValue + mu * (t - currentTime), sigma * (t - currentTime));
}

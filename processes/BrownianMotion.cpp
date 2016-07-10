#include "BrownianMotion.h"

BrownianMotion::BrownianMotion(double drift, double volatility, double deltaT) :
    StableProcess(2.0, 0.0, drift, volatility * M_SQRT1_2, deltaT)
{
}

void BrownianMotion::nextImpl()
{
    currentValue += mu * dt + sigma * NormalRand::variate(0, dtCoef);
}

void BrownianMotion::nextImpl(double deltaT)
{
    currentValue += mu * deltaT + sigma * NormalRand::variate(0, std::sqrt(deltaT));
}

double BrownianMotion::QuantileImpl(double t, double p) const
{
    return NormalRand::quantile(p, currentValue, std::sqrt(t - currentTime));
}

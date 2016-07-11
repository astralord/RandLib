#include "BrownianMotion.h"

BrownianMotion::BrownianMotion(double deltaT) :
    StableProcess(2.0, 0.0, 0.0, M_SQRT1_2, deltaT)
{

}

BrownianMotion::BrownianMotion(double drift, double volatility, double deltaT) :
    StableProcess(2.0, 0.0, drift, volatility * M_SQRT1_2, deltaT)
{
}

void BrownianMotion::nextImpl()
{
    currentValue += mu * dt + sigma * NormalRand::variate(0, dtCoef);
}

double BrownianMotion::QuantileImpl(double t, double p) const
{
    return NormalRand::quantile(p, currentValue + mu * (t - currentTime), sigma * std::sqrt(t - currentTime));
}

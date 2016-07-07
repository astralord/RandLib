#include "BrownianMotion.h"

BrownianMotion::BrownianMotion(double deltaT) :
    StochasticProcess(deltaT),
    sqrtDt(std::sqrt(dt))
{
    currentValue = 0;
}

double BrownianMotion::nextImpl()
{
    currentValue += NormalRand::variate(0, sqrtDt);
    return currentValue;
}

double BrownianMotion::nextImpl(double deltaT)
{
    currentValue += NormalRand::variate(0, std::sqrt(deltaT));
    return currentValue;
}

double BrownianMotion::MeanImpl(double) const
{
    return currentValue;
}

double BrownianMotion::VarianceImpl(double t) const
{
    return t - currentTime;
}

double BrownianMotion::QuantileImpl(double t, double p) const
{
    return NormalRand::quantile(p, currentValue, std::sqrt(t - currentTime));
}

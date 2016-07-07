#include "PoissonProcess.h"

PoissonProcess::PoissonProcess(double deltaT, double rate) :
    StochasticProcess(deltaT),
    lambda(rate > 0.0 ? rate : 1.0),
    futureJumpTime(ExponentialRand::variate(lambda))
{

}

void PoissonProcess::nextImpl()
{
    while (currentTime > futureJumpTime) {
        ++currentValue;
        futureJumpTime += ExponentialRand::variate(lambda);
    }
}

void PoissonProcess::nextImpl(double)
{
    nextImpl();
}

double PoissonProcess::MeanImpl(double t) const
{
    return currentValue + lambda * (t - currentTime);
}

double PoissonProcess::VarianceImpl(double t) const
{
    return lambda * (t - currentTime);
}

double PoissonProcess::QuantileImpl(double t, double p) const
{
    PoissonRand X(lambda * (t - currentTime));
    return X.Quantile(p) + currentValue;
}


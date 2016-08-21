#include "PoissonProcess.h"

PoissonProcess::PoissonProcess(double rate, double deltaT) :
    StochasticProcess(deltaT),
    lambda(rate > 0.0 ? rate : 1.0),
    futureJumpTime(ExponentialRand::Variate(lambda))
{

}

void PoissonProcess::nextImpl()
{
    while (currentTime > futureJumpTime) {
        ++currentValue;
        futureJumpTime += ExponentialRand::Variate(lambda);
    }
}

double PoissonProcess::MeanImpl(double t) const
{
    return currentValue + lambda * (t - currentTime);
}

double PoissonProcess::VarianceImpl(double t) const
{
    return lambda * (t - currentTime);
}

double PoissonProcess::Quantile(double t, double p) const
{
    PoissonRand X(lambda * (t - currentTime));
    return X.Quantile(p) + currentValue;
}


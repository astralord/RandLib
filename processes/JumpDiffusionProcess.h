#ifndef JUMPDIFFUSIONPROCESS_H
#define JUMPDIFFUSIONPROCESS_H

#include "StochasticProcess.h"
#include "BrownianMotion.h"
#include "CompoundPoissonProcess.h"

/**
 * @brief The JumpDiffusionProcess class
 * dX(t) = mu * dt + sigma * dB(t) + dJ(t),
 * where B(t) is a Brownian motion, J(t) is a Compound poisson process,
 * mu is drift and sigma is volatility.
 */
template <typename T>
class RANDLIBSHARED_EXPORT JumpDiffusionProcess : public StochasticProcess
{
    BrownianMotion B;
    CompoundPoissonProcess<T> J;

public:
    JumpDiffusionProcess(double drift, double volatility, double rate, const UnivariateProbabilityDistribution<T> &jumpDistribution, double deltaT = 1.0);

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};

#endif // JUMPDIFFUSIONPROCESS_H

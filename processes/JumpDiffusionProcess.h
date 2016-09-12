#ifndef JUMPDIFFUSIONPROCESS_H
#define JUMPDIFFUSIONPROCESS_H

#include "BrownianMotion.h"
#include "CompoundPoissonProcess.h"

/**
 * @brief The JumpDiffusionProcess class
 * dX(t) = μdt + σdB(t) + dJ(t),
 * where B(t) is a Brownian motion, J(t) is a Compound Poisson process.
 * μ is drift and σ is volatility.
 *
 * Notation: JD(t | μ, σ, λ, F)
 *
 * Closed-form solution:
 * JD(t | μ, σ, λ, F) = B(t | μ, σ) + J(t | λ, F)
 */
template <typename T>
class RANDLIBSHARED_EXPORT JumpDiffusionProcess : public StochasticProcess<double>
{
    BrownianMotion B;
    CompoundPoissonProcess<T> J;

public:
    JumpDiffusionProcess(double drift, double volatility, double rate, const UnivariateProbabilityDistribution<T> &jumpDistribution, double deltaT = 1.0);

private:
    void NextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
};

#endif // JUMPDIFFUSIONPROCESS_H

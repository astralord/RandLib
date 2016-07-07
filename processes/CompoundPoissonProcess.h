#ifndef COMPOUNDPOISSONPROCESS_H
#define COMPOUNDPOISSONPROCESS_H

#include "PoissonProcess.h"
#include "../distributions/univariate/UnivariateProbabilityDistribution.h"

/**
 * @brief The CompoundPoissonProcess class
 */
template <typename T>
class RANDLIBSHARED_EXPORT CompoundPoissonProcess : public StochasticProcess
{
    PoissonProcess N;
    const UnivariateProbabilityDistribution<T> &Y;
    int jumpsAmount;
public:
    CompoundPoissonProcess(double rate, const UnivariateProbabilityDistribution<T> &jumpDistribution, double deltaT = 1.0);

private:
    void nextImpl() override;
    void nextImpl(double deltaT) override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};

#endif // COMPOUNDPOISSONPROCESS_H

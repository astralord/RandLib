#ifndef POISSONPROCESS_H
#define POISSONPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/discrete/PoissonRand.h"
#include "../distributions/univariate/continuous/ExponentialRand.h"

/**
 * @brief The PoissonProcess class
 * X(t) = #(S_j < t)
 * where S_j = \sum_1^j W_i,
 * where W_i ~ Exp(λ)
 *
 * Notation: N(t | λ)
 */
class RANDLIBSHARED_EXPORT PoissonProcess : public StochasticProcess<int>
{
    double lambda;
    double futureJumpTime;
public:
    explicit PoissonProcess(double rate, double deltaT = 1.0);

    inline double GetRate() const { return lambda; }

private:
    void NextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
public:
    double Quantile(double t, double p) const;
};

#endif // POISSONPROCESS_H

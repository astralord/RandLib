#ifndef POISSONPROCESS_H
#define POISSONPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/discrete/PoissonRand.h"
#include "../distributions/univariate/continuous/ExponentialRand.h"

class RANDLIBSHARED_EXPORT PoissonProcess : public StochasticProcess
{
    double lambda;
    double futureJumpTime;
public:
    explicit PoissonProcess(double deltaT = 1.0, double rate = 1.0);

private:
    void nextImpl() override;
    void nextImpl(double deltaT) override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};

#endif // POISSONPROCESS_H

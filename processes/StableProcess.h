#ifndef STABLEPROCESS_H
#define STABLEPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/continuous/StableRand.h"

class RANDLIBSHARED_EXPORT StableProcess : public StochasticProcess
{
    StableRand X;
    double dtCoef;
    double drift, volatility;

public:
    explicit StableProcess(double exponent, double skewness, double scale, double location, double deltaT = 1.0);

private:
    void nextImpl() override;
    void nextImpl(double deltaT) override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};

#endif // STABLEPROCESS_H

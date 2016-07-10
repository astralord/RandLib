#ifndef STABLEPROCESS_H
#define STABLEPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/continuous/StableRand.h"

class RANDLIBSHARED_EXPORT StableProcess : public StochasticProcess
{
protected:
    StableRand X;
    double dtCoef;
    double mu, sigma;

public:
    StableProcess(double exponent, double skewness, double drift, double volatility, double deltaT = 1.0);

private:
    void nextImpl() override;
    void nextImpl(double deltaT) override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
    void QuantileImpl(const std::vector<double> &t, std::vector<double> &outputData, double p) const;
};

#endif // STABLEPROCESS_H

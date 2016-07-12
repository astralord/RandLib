#ifndef STABLEPROCESS_H
#define STABLEPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/continuous/StableRand.h"

/**
 * @brief The StableProcess class
 * dX(t) = μdt + σdS(t),
 * where S(t) - S(s) ~ Stable(α, β, (t - s)^(1/α), 0),
 * α is characteristic exponent, β is skewness, μ is drift and σ is volatility.
 * α == 2 - X(t) is a Brownian motion (Wiener process)
 * α == 1, beta == 0 - X(t) is a Cauchy process
 */
class RANDLIBSHARED_EXPORT StableProcess : public StochasticProcess<double>
{
protected:
    StableRand X;
    double dtCoef;
    double mu, sigma;

public:
    StableProcess(double exponent, double skewness, double drift, double volatility, double deltaT = 1.0);

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
    void QuantileImpl(const std::vector<double> &t, std::vector<double> &outputData, double p) const;
};

#endif // STABLEPROCESS_H

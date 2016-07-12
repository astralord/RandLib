#ifndef ORNSTEIN_UHLENBECKPROCESS_H
#define ORNSTEIN_UHLENBECKPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/continuous/NormalRand.h"

/**
 * @brief The OrnsteinUhlenbeckProcess class
 * dX(t) = (alpha - beta * X(t)) * dt + sigma * dB(t),
 * where B(t) is a Brownian motion, alpha is drift, beta is reversion speed and sigma is volatility.
 */
class RANDLIBSHARED_EXPORT OrnsteinUhlenbeckProcess : public StochasticProcess<double>
{
    double alpha, beta, sigma;
    double expmBetaDt;
    NormalRand X;

public:
    OrnsteinUhlenbeckProcess(double drift, double reversionSpeed, double volatility, double initialValue, double deltaT = 1.0);

    inline double getDrift() { return alpha; }
    inline double getReversionSpeed() { return beta; }
    inline double getVolatility() { return sigma; }

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};

#endif // ORNSTEIN_UHLENBECKPROCESS_H
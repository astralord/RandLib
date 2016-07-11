#ifndef GEOMETRICBROWNIANMOTION_H
#define GEOMETRICBROWNIANMOTION_H

#include "StochasticProcess.h"
#include "BrownianMotion.h"
#include "../distributions/univariate/continuous/LogNormalRand.h"

/**
 * @brief The GeometricBrownianMotion class
 * dX(t) = X(t) * (mu * dt + sigma * dB(t)),
 * where B(t) is a Brownian motion, mu is drift and sigma is volatility.
 */
class RANDLIBSHARED_EXPORT GeometricBrownianMotion : public StochasticProcess
{
    double mu, sigma, S0;
    double mumSigma2_2;
    BrownianMotion B;

public:
    GeometricBrownianMotion(double drift, double volatility, double initialValue, double deltaT = 1.0);

    inline double getDrift() { return mu; }
    inline double getVolatility() { return sigma; }
    inline double getInitialValue() { return S0; }

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;

    double QuantileImpl(double t, double p) const override;
};

#endif // GEOMETRICBROWNIANMOTION_H

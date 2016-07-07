#ifndef GEOMETRICBROWNIANMOTION_H
#define GEOMETRICBROWNIANMOTION_H

#include "StochasticProcess.h"
#include "BrownianMotion.h"
#include "../distributions/univariate/continuous/LogNormalRand.h"

class RANDLIBSHARED_EXPORT GeometricBrownianMotion : public StochasticProcess
{
    double mu, sigma, S0;
    double mumSigma2_2;
    BrownianMotion B;

public:
    GeometricBrownianMotion(double deltaT = 1.0, double drift = 0.0, double volatility = 1.0, double initialValue = 1.0);

    inline double getDrift() { return mu; }
    inline double getVolatility() { return sigma; }
    inline double getInitialValue() { return S0; }

private:
    double nextImpl() override;
    double nextImpl(double deltaT) override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;

    double QuantileImpl(double t, double p) const override;
};

#endif // GEOMETRICBROWNIANMOTION_H

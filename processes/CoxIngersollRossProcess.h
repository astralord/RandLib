#ifndef COXINGERSOLLROSSPROCESS_H
#define COXINGERSOLLROSSPROCESS_H

#include "BrownianMotion.h"
#include "../distributions/univariate/continuous/NoncentralChiSquared.h"

/**
 * @brief The CoxIngersollRossProcess class
 * dX(t) = (α - βX(t))dt + σ√(X(t))dB(t),
 * where B(t) is a Brownian motion, α is drift, β is reversion speed and σ is volatility.
 */
class RANDLIBSHARED_EXPORT CoxIngersollRossProcess : public StochasticProcess<double>
{
    double alpha, beta, sigma;
    double expmBetaDt, c, degree;

public:
    CoxIngersollRossProcess(double drift, double reversionSpeed, double volatility, double initialValue, double deltaT = 1.0);

    inline double getDrift() { return alpha; }
    inline double getReversionSpeed() { return beta; }
    inline double getVolatility() { return sigma; }

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};


#endif // COXINGERSOLLROSSPROCESS_H

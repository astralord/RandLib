#ifndef HESTONPROCESS_H
#define HESTONPROCESS_H

#include "CoxIngersollRossProcess.h"
#include "../distributions/univariate/continuous/LogNormalRand.h"

/**
 * @brief The HestonProcess class
 * dX(t) = μX(t)dt + √(V(t))X(t)dW(t)),
 * where V(t) is CIR process (volatility) satisfying:
 * dV(t) = (α - β V(t)) dt + σ √(V(t)) dB(t).
 * B(t) and W(t) are correlated Wiener processes with dW(t) dB(t) = ρ dt,
 * μ is the drift term (the rate of return of the asset),
 * α is the drift of V(t),
 * β is the reversion speed of V(t),
 * σ is the volatility of V(t),
 * ρ is the correlation between B(t) and W(t).
 */
class RANDLIBSHARED_EXPORT HestonProcess : public StochasticProcess<double>
{
    double mu, alpha, beta, sigma, ro;
    CoxIngersollRossProcess V;

public:
    HestonProcess(double drift, double volatilityDrift, double reversionSpeed, double volatility, double initialValue, double volatilityInitialValue, double correlation, double deltaT = 1.0);

    inline double getDrift() { return mu; }
    inline double getVolatilityDrift() { return alpha; }
    inline double getReversionSpeed() { return beta; }
    inline double getVolatility() { return sigma; }
    inline double getCorrelation() { return ro; }

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};

#endif // HESTONPROCESS_H

#ifndef COXINGERSOLLROSSPROCESS_H
#define COXINGERSOLLROSSPROCESS_H

#include "BrownianMotion.h"
#include "../distributions/univariate/continuous/NoncentralChiSquared.h"

/**
 * @brief The CoxIngersollRossProcess class
 * dX(t) = β(α - X(t))dt + σ(X(t))^(1/2)dB(t),
 * where B(t) is a Brownian motion, α is mean, β is reversion speed and σ is volatility.
 *
 * Notation: CIR(t | α, β, σ, X0)
 *
 * Closed-form solution:
 * CIR(t | α, β, σ, X0) ~ Y / c, where
 * Y ~ Noncentral Chi-Squared(4αβ / σ^2, X0 * c * exp(-βt)))
 * and c = 4β / (σ^2 * (1.0 - exp(-βt))).
 *
 * Special cases:
 * CIR(t | α, β, 0, X0) = α + (X0 - α) exp(-βt)
 */
class RANDLIBSHARED_EXPORT CoxIngersollRossProcess : public StochasticProcess<double>
{
    double alpha, beta, sigma;
    double expmBetaDt, c, degree;

public:
    CoxIngersollRossProcess(double mean, double reversionSpeed, double volatility, double initialValue, double deltaT = 1.0);

    inline double GetMean() { return alpha; }
    inline double GetReversionSpeed() { return beta; }
    inline double GetVolatility() { return sigma; }

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
};


#endif // COXINGERSOLLROSSPROCESS_H

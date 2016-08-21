#ifndef ORNSTEIN_UHLENBECKPROCESS_H
#define ORNSTEIN_UHLENBECKPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/continuous/NormalRand.h"

/**
 * @brief The OrnsteinUhlenbeckProcess class
 * dX(t) = β(α - X(t))dt + σdB(t),
 * where B(t) is a standard Brownian motion,
 * α is mean, β is reversion speed and σ is volatility.
 *
 * Notation: OU(t | α, β, σ, X0)
 *
 * Closed-form solution:
 * OU(t | α, β, σ) = X0 exp(-βt) + Y,
 * where Y ~ Normal(α(1 - exp(-βt)), σ^2 (1 - exp(-2βt)) / 2β)
 *
 * Special cases:
 * OU(t | α, β, 0) = X0 exp(-βt) + α(1 - exp(-βt)),
 * OU(t | α, 0, σ) = X0 + B(t | 0, σ)
 */
class RANDLIBSHARED_EXPORT OrnsteinUhlenbeckProcess : public StochasticProcess<double>
{
    double alpha, beta, sigma;
    double expmBetaDt;
    NormalRand X;

public:
    OrnsteinUhlenbeckProcess(double mean, double reversionSpeed, double volatility, double initialValue, double deltaT = 1.0);

    inline double GetMean() { return alpha; }
    inline double GetReversionSpeed() { return beta; }
    inline double GetVolatility() { return sigma; }

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;

public:
    double Quantile(double t, double p) const;
};

#endif // ORNSTEIN_UHLENBECKPROCESS_H

#ifndef STABLEPROCESS_H
#define STABLEPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/continuous/StableRand.h"

/**
 * @brief The StableProcess class
 * dX(t) = μdt + σdS(t),
 * where S(t) - S(s) ~ Stable(α, β, (t - s)^(1/α), 0),
 * α is characteristic exponent, β is skewness, μ is drift and σ is volatility.
 *
 * Notation: Stable(t | α, β, μ, σ).
 *
 * Closed-form solution:
 * Stable(t | α, β, μ, σ) = μt + σ t^(1/α) S,
 * where S ~ Stable(α, β, 1, 0)
 *
 * Special cases:
 * Stable(t | 2, β, μ, σ) = B(t | μ, σ),
 * Stable(t | 1, 0, μ, σ) = Cauchy(t | μ, σ).
 */
class RANDLIBSHARED_EXPORT StableProcess : public StochasticProcess<double>
{
protected:
    StableRand X;
    double dtCoef;
    double mu, sigma;

public:
    StableProcess(double exponent, double skewness, double drift = 0.0, double volatility = 1.0, double deltaT = 1.0);

    double Quantile(double t, double p) const;
    void Quantile(const std::vector<double> &t, std::vector<double> &outputData, double p) const;
private:
    void NextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
};

#endif // STABLEPROCESS_H

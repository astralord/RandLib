#ifndef GEOMETRICBROWNIANMOTION_H
#define GEOMETRICBROWNIANMOTION_H

#include "BrownianMotion.h"
#include "../distributions/univariate/continuous/LogNormalRand.h"

/**
 * @brief The GeometricBrownianMotion class
 * dX(t) = μX(t)dt + σX(t)dB(t),
 * where B(t) is a Brownian motion, μ is drift and σ is volatility.
 *
 * Notation: GBM(t | μ, σ).
 *
 * Closed-form solution:
 * GBM(t | μ, σ) = exp((μ - σ^2/2)t - σB(t))
 *
 * Special cases:
 * GBM(t | μ, 0) = exp(μt).
 */
class RANDLIBSHARED_EXPORT GeometricBrownianMotion : public StochasticProcess<double>
{
    double mu, sigma;
    double mumSigma2_2;
    LogNormalRand X;

public:
    GeometricBrownianMotion(double drift, double volatility, double initialValue, double deltaT = 1.0);

    inline double getDrift() { return mu; }
    inline double getVolatility() { return sigma; }

private:
    void nextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;

public:
    void ProbabilityDensityFunction(double t, const std::vector<double> &x, std::vector<double> &y) const;
    void CumulativeDistributionFunction(double t, const std::vector<double> &x, std::vector<double> &y) const;
    double Quantile(double t, double p) const;
};

#endif // GEOMETRICBROWNIANMOTION_H

#ifndef BROWNIANMOTION_H
#define BROWNIANMOTION_H

#include "StableProcess.h"
#include "../distributions/univariate/continuous/NormalRand.h"

/**
 * @brief The BrownianMotion class
 * dX(t) = μdt + σdB(t),
 * where B(t) - B(s) ~ Normal(0, t - s) for t > s, μ is drift and σ is volatility.
 *
 * Notation: B(t | μ, σ).
 */
class RANDLIBSHARED_EXPORT BrownianMotion : public StableProcess
{
public:
    explicit BrownianMotion(double deltaT = 1.0);
    BrownianMotion(double drift, double volatility, double deltaT = 1.0);

private:
    void NextImpl() override;

public:
    void ProbabilityDensityFunction(double t, const std::vector<double> &x, std::vector<double> &y) const;
    void CumulativeDistributionFunction(double t, const std::vector<double> &x, std::vector<double> &y) const;
    double Quantile(double t, double p) const;
};


#endif // BROWNIANMOTION_H

#ifndef CAUCHYPROCESS_H
#define CAUCHYPROCESS_H

#include "StableProcess.h"

/**
 * @brief The CauchyProcess class
 * dX(t) = μdt + σdC(t),
 * where C(t) - C(s) ~ Cauchy(0, t - s) for t > s, μ is drift and σ is volatility.
 *
 * Notation: Cauchy(t | μ, σ).
 */
class RANDLIBSHARED_EXPORT CauchyProcess : public StableProcess
{
public:
    CauchyProcess(double drift = 0.0, double volatility = 1.0, double deltaT = 1.0);

private:
    void NextImpl() override;

public:
    void ProbabilityDensityFunction(double t, const std::vector<double> &x, std::vector<double> &y) const;
    void CumulativeDistributionFunction(double t, const std::vector<double> &x, std::vector<double> &y) const;
    double Quantile(double t, double p) const;
};


#endif // CAUCHYPROCESS_H

#ifndef CAUCHYPROCESS_H
#define CAUCHYPROCESS_H

#include "StableProcess.h"

/**
 * @brief The CauchyProcess class
 * dX(t) = μdt + σdC(t),
 * where C(t) - C(s) ~ Cauchy(0, t - s) for t > s, μ is drift and σ is volatility.
 */
class RANDLIBSHARED_EXPORT CauchyProcess : public StableProcess
{
public:
    CauchyProcess(double drift, double volatility, double deltaT = 1.0);
    void nextImpl() override;

    double QuantileImpl(double t, double p) const override;
};


#endif // CAUCHYPROCESS_H

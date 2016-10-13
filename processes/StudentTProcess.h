#ifndef STUDENTTPROCESS_H
#define STUDENTTPROCESS_H

#include "StochasticProcess.h"
#include "../distributions/univariate/continuous/StudentTRand.h"

/**
 * @brief The StudentTProcess class
 * dX(t) = μdt + σdB(t),
 * where B(t) - B(s) ~ Normal(0, t - s) for t > s, μ is drift and σ is volatility.
 *
 * Notation: B(t | μ, σ).
 */
class RANDLIBSHARED_EXPORT StudentTProcess : public StochasticProcess<double>
{
    StudentTRand X;
    double mu, sigma;

public:
    StudentTProcess(int degree, double drift = 0, double volatility = 1.0, double deltaT = 1.0);

    inline double GetDegree() { return X.GetDegree(); }
    inline double GetDrift() { return mu; }
    inline double GetVolatility() { return sigma; }

private:
    void NextImpl() override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
};


#endif // STUDENTTPROCESS_H

#ifndef STABLERAND_H
#define STABLERAND_H

#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"
#include "LevyRand.h"

class RANDLIBSHARED_EXPORT StableRand : public ContinuousRand
{
    double alpha, beta, sigma, mu;
    UniformRand U;
    ExponentialRand E;
    double B, S, alphaInv; /// coefficients for alpha != 1
    double logSigma; /// coefficients for alpha == 1

public:
    StableRand(double exponent, double skewness, double scale, double location);

    void setAlphaAndBeta(double exponent, double skewness);
    void setSigma(double scale);
    void setMu(double location);

    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }
    inline double getSigma()  const{ return sigma; }
    inline double getMu() const { return mu; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    inline double M() const override { return (alpha > 1) ? mu : NAN; }
    inline double Var() const override { return (alpha == 2) ? 2 * sigma * sigma : INFINITY; }
};



#endif // STABLERAND_H

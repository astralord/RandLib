#ifndef STABLERAND_H
#define STABLERAND_H

#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"
#include "LevyRand.h"

class RANDLIBSHARED_EXPORT StableRand : public ContinuousRandom
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

    inline double getAlpha() { return alpha; }
    inline double getBeta() { return beta; }
    inline double getSigma() { return sigma; }
    inline double getMu() { return mu; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return (alpha > 1) ? mu : NAN; }
    inline double Var() { return (alpha == 2) ? 2 * sigma * sigma : INFINITY; }
};



#endif // STABLERAND_H

#ifndef BETARAND_H
#define BETARAND_H

#include <RandomVariable.h>
#include "GammaRand.h"

class RANDLIBSHARED_EXPORT BetaRand : public ContinuousRandom
{
    double alpha, beta;
    GammaRand X, Y;
    double gammaA, gammaB, pdfCoef;
public:
    BetaRand(double shape1, double shape2);

    void setParameters(double shape1, double shape2);
    void setAlpha(double shape1);
    void setBeta(double shape2);
    double getAlpha() { return alpha; }
    double getBeta() { return beta; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return alpha / (alpha + beta); }
    inline double Var() {
        double denom = alpha + beta;
        denom *= denom * (denom + 1);
        return alpha * beta / denom;
    }
};

#endif // BETARAND_H

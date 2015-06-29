#ifndef BETARAND_H
#define BETARAND_H

#include "ContinuousRand.h"
#include "GammaRand.h"

class RANDLIBSHARED_EXPORT BetaRand : public ContinuousRand
{
    double alpha, beta;
    GammaRand X, Y;
    double gammaA, gammaB, pdfCoef;
public:
    BetaRand(double shape1, double shape2);

    void setParameters(double shape1, double shape2);
    void setAlpha(double shape1);
    void setBeta(double shape2);
    inline double getAlpha() const { return alpha; }
    inline double getBeta() const { return beta; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double M() const override { return alpha / (alpha + beta); }
    double Var() const override {
        double denom = alpha + beta;
        denom *= denom * (denom + 1);
        return alpha * beta / denom;
    }
};

#endif // BETARAND_H

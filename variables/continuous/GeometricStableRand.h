#ifndef GEOMETRICSTABLERAND_H
#define GEOMETRICSTABLERAND_H

#include "StableRand.h"
#include "ExponentialRand.h"

class RANDLIBSHARED_EXPORT GeometricStableRand : public ContinuousRand
{
    double mu, sigma;

    StableRand S;

public:
    GeometricStableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual void setName() override;

    void setParameters(double exponent, double skewness, double scale, double location);
    inline double getAlpha() const { return S.getAlpha(); }
    inline double getBeta() const { return S.getBeta(); }
    inline double getSigma() const { return sigma; }
    inline double getMu() const { return mu; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return (S.getAlpha()> 1) ? mu : NAN; }
    // write correct Variance!
    double Var() const override { return (S.getAlpha() == 2) ? sigma : INFINITY; }
};

#endif // GEOMETRICSTABLERAND_H

#ifndef GAMMARAND_H
#define GAMMARAND_H

#include <RandomVariable.h>
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

class RANDLIBSHARED_EXPORT GammaRand : public ContinuousRand
{
    double k, theta;
    double kInv, thetaInv; /// 1.0 / k and 1.0 / theta
    double valueCoef; /// (e + k) / (k * e)
    double cdfCoef; /// 1.0 / gamma(k)
    double pdfCoef; /// 1.0 / (gamma(k) * theta ^ k)
    UniformRand U;
    ExponentialRand E;
    NormalRand N;

    double m, s, s_2, d, b, w, v, c; /// constants for sampling
    void setConstants();

public:
    GammaRand(double shape, double scale);

    void setParameters(double shape, double scale);
    void setShape(double shape);
    void setScale(double scale);
    inline double getShape() const { return k; }
    inline double getScale() const { return theta; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    inline double M() const override { return k * theta; }
    inline double Var() const override { return k * theta * theta; }
};

#endif // GAMMARAND_H

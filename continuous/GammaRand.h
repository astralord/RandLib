#ifndef GAMMARAND_H
#define GAMMARAND_H

#include <RandomVariable.h>

class RANDLIBSHARED_EXPORT GammaRand : public ContinuousRand
{
    double k, theta;
    double thetaInv; // 1.0 / theta
    double cdfCoef; // 1.0 / gamma(k)
    double pdfCoef; // 1.0 / (gamma(k) * theta ^ k)
public:
    GammaRand(double shape, double scale);

    void setParameters(double shape, double scale);
    void setShape(double shape);
    void setScale(double scale);
    double getShape() { return k; }
    double getScale() { return theta; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return k * theta; }
    inline double Var() { return k * theta * theta; }
};

#endif // GAMMARAND_H

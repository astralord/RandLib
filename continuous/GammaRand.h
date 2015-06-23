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
    inline double getShape() const { return k; }
    inline double getScale() const { return theta; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    inline double M() const override { return k * theta; }
    inline double Var() const override { return k * theta * theta; }
};

#endif // GAMMARAND_H

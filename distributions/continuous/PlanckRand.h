#ifndef PLANCKRAND_H
#define PLANCKRAND_H

#include "GammaRand.h"
#include "../discrete/ZetaRand.h"

class RANDLIBSHARED_EXPORT PlanckRand : public ContinuousDistribution
{
    double a, b;
    double pdfCoef;

    ZetaRand Z;
    GammaRand G;

public:
    PlanckRand(double shape, double scale);
    std::string name() override;

    void setParameters(double shape, double scale);
    inline double getShape() const { return a; }
    inline double getScale() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;
};

#endif // PLANCKRAND_H

#ifndef PLANCKRAND_H
#define PLANCKRAND_H

#include "GammaRand.h"
#include "../discrete/ZetaRand.h"

/**
 * @brief The PlanckRand class
 */
class RANDLIBSHARED_EXPORT PlanckRand : public ContinuousDistribution
{
    double a, b;
    double pdfCoef;

    ZetaRand Z;
    GammaRand G;

public:
    PlanckRand(double shape, double scale);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void setParameters(double shape, double scale);
    inline double getShape() const { return a; }
    inline double getScale() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
};

#endif // PLANCKRAND_H

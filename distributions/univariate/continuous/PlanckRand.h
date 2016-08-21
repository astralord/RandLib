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

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double scale);
    inline double GetShape() const { return a; }
    inline double GetScale() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
};

#endif // PLANCKRAND_H

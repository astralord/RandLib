#ifndef INVERSEGAMMARAND_H
#define INVERSEGAMMARAND_H

#include "GammaRand.h"

/**
 * @brief The InverseGammaRand class
 * Inverse-Gamma distribution
 * X ~ Inv-Gamma(\alpha, \beta)
 */
class RANDLIBSHARED_EXPORT InverseGammaRand : public GammaRand
{
public:
    InverseGammaRand(double shape = 1, double scale = 1);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // INVERSEGAMMARAND_H

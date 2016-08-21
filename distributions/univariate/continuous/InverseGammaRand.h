#ifndef INVERSEGAMMARAND_H
#define INVERSEGAMMARAND_H

#include "GammaRand.h"

/**
 * @brief The InverseGammaRand class
 * Inverse-Gamma distribution
 * X ~ Inv-Gamma(α, β)
 *
 * X = 1 / Y, where Y ~ Gamma(α, β)
 */
class RANDLIBSHARED_EXPORT InverseGammaRand : public ContinuousDistribution
{
    GammaRand X;
    double alpha, beta;
public:
    InverseGammaRand(double shape = 1, double rate = 1);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }
    void setParameters(double shape, double rate);
    inline double getShape() const { return alpha; }
    inline double getRate() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    double getLogGammaFunction() const { return X.getLogGammaFunction(); }
};

#endif // INVERSEGAMMARAND_H

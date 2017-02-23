#ifndef INVERSEGAMMARAND_H
#define INVERSEGAMMARAND_H

#include "GammaRand.h"

/**
 * @brief The InverseGammaRand class
 * Inverse-Gamma distribution
 *
 * X ~ Inv-Gamma(α, β)
 *
 * Related distributions:
 * X = 1 / Y, where Y ~ Gamma(α, β)
 */
class RANDLIBSHARED_EXPORT InverseGammaRand : public ContinuousDistribution
{
    GammaRand X;
    double alpha, beta;
    double pdfCoef;
public:
    InverseGammaRand(double shape = 1, double rate = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }
    void SetParameters(double shape, double rate);
    inline double GetShape() const { return alpha; }
    inline double GetRate() const { return beta; }

    double f(double x) const override;
    double logf(double x) const override;
    double F(double x) const override;
    double S(double x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    double GetLogGammaFunction() const { return X.GetLogGammaFunction(); }
};

#endif // INVERSEGAMMARAND_H

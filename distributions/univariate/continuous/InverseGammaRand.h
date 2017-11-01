#ifndef INVERSEGAMMARAND_H
#define INVERSEGAMMARAND_H

#include "GammaRand.h"

/**
 * @brief The InverseGammaRand class <BR>
 * Inverse-Gamma distribution
 *
 * X ~ Inv-Gamma(α, β)
 *
 * Related distributions: <BR>
 * X = 1 / Y, where Y ~ Gamma(α, β)
 */
class RANDLIBSHARED_EXPORT InverseGammaRand : public ContinuousDistribution
{
    double alpha = 1; ///< shape α
    double beta = 1; ///< rate β
    double pdfCoef = 0; ///< coefficient for faster pdf calculation

    GammaRand X{};
public:
    InverseGammaRand(double shape = 1, double rate = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }
    void SetParameters(double shape, double rate);
    inline double GetShape() const { return alpha; }
    inline double GetRate() const { return beta; }
    inline double GetLogShape() const { return X.GetLogShape(); }
    inline double GetLogRate() const { return X.GetLogRate(); }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
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
    double GetLogGammaShape() const { return X.GetLogGammaShape(); }
};

#endif // INVERSEGAMMARAND_H

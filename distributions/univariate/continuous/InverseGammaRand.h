#ifndef INVERSEGAMMARAND_H
#define INVERSEGAMMARAND_H

#include "GammaRand.h"

/**
 * @brief The InverseGammaRand class <BR>
 * Inverse-Gamma distribution
 *
 * X ~ Inv-Γ(α, β)
 *
 * Related distributions: <BR>
 * X = 1 / Y, where Y ~ Gamma(α, β)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT InverseGammaRand : public ContinuousDistribution<RealType>
{
    double alpha = 1; ///< shape α
    double beta = 1; ///< rate β
    double pdfCoef = 0; ///< coefficient for faster pdf calculation

    GammaRand<RealType> X{};

public:
    InverseGammaRand(double shape = 1, double rate = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double rate);
    inline double GetShape() const { return alpha; }
    inline double GetRate() const { return beta; }
    inline double GetLogShape() const { return X.GetLogShape(); }
    inline double GetLogRate() const { return X.GetLogRate(); }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

public:
    double GetLogGammaShape() const { return X.GetLogGammaShape(); }
};

#endif // INVERSEGAMMARAND_H

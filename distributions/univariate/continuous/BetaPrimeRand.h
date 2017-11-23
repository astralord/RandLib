#ifndef BETAPRIMERAND_H
#define BETAPRIMERAND_H

#include "BetaRand.h"

/**
 * @brief The BetaPrimeRand class <BR>
 * Beta-prime distribution
 *
 * f(x | α, β) = x^{α-1} (1 + x)^{-α - β} / B(α, β), <BR>
 * where B(α, β) denotes Beta function
 *
 * Notation: X ~ B'(α, β)
 *
 * Related distributions: <BR>
 * X / (X + 1) ~ B(α, β) <BR>
 * X = Y / Z, where Y ~ Γ(α) and Z ~ Γ(β) <BR>
 * β/α * X ~ F(2α, 2β)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT BetaPrimeRand : public ContinuousDistribution<RealType>
{
    double alpha = 1; ///< first shape α
    double beta = 1; ///< second shape β
    BetaRand<RealType> B{};

public:
    BetaPrimeRand(double shape1 = 1, double shape2 = 1);
    String Name() const override;
    void SetShapes(double shape1, double shape2);
    inline double GetAlpha() const { return alpha; }
    inline double GetBeta() const { return beta; }

    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

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

    std::complex<double> CFImpl(double t) const override;
public:
    /**
     * @fn GetBetaFunction
     * @return B(α, β)
     */
    inline double GetBetaFunction() const { return B.GetBetaFunction(); }

    /**
     * @fn GetLogBetaFunction
     * @return log(B(α, β))
     */
    inline double GetLogBetaFunction() const { return B.GetLogBetaFunction(); }

    /**
     * @fn FitAlpha
     * fit α by maximum-likelihood
     * @param sample
     */
    void FitAlpha(const std::vector<RealType> &sample);

    /**
     * @fn FitBeta
     * fit β by maximum-likelihood
     * @param sample
     */
    void FitBeta(const std::vector<RealType> &sample);

    /**
     * @fn Fit
     * fit shapes by maximum-likelihood
     * @param sample
     */
    void Fit(const std::vector<RealType> &sample);
};

#endif // BETAPRIMERAND_H

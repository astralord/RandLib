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
 * X ~ F(2α, 2β)
 */
class RANDLIBSHARED_EXPORT BetaPrimeRand : public ContinuousDistribution
{
    double alpha = 1; ///< first shape α
    double beta = 1; ///< second shape β
    BetaRand B{};

public:
    BetaPrimeRand(double shape1 = 1, double shape2 = 1);
    String Name() const override;
    void SetParameters(double shape1, double shape2);
    inline double GetAlpha() const { return alpha; }
    inline double GetBeta() const { return beta; }

    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;

    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

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
};

#endif // BETAPRIMERAND_H

#ifndef BETAPRIMERAND_H
#define BETAPRIMERAND_H

#include "BetaRand.h"

/**
 * @brief The BetaPrimeRand class
 * Beta-prime distribution
 *
 * f(x | α, β) = x^{α-1} (1 + x)^{-α - β} / B(α, β),
 * where B(α, β) denotes Beta function
 *
 * Notation: X ~ Beta-Prime(α, β)
 *
 * Related distributions:
 * X / (X + 1) ~ Beta(α, β)
 * X = Y / Z, where Y ~ Gamma(α) and Z ~ Gamma(β)
 * X ~ F(2α, 2β)
 */
class RANDLIBSHARED_EXPORT BetaPrimeRand : public ContinuousDistribution
{
    double alpha, beta;
    BetaRand B;
public:
    BetaPrimeRand(double shape1 = 1, double shape2 = 1);
    std::string Name() const override;
    void SetParameters(double shape1, double shape2);
    inline double GetAlpha() const { return alpha; }
    inline double GetBeta() const { return beta; }

    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

    void Sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CFImpl(double t) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
public:
    /**
     * @brief GetInvBetaFunction
     * @return 1 / Beta(α, β)
     */
    inline double GetInverseBetaFunction() const { return B.GetInverseBetaFunction(); }

    /**
     * @brief GetLogBetaFunction
     * @return log Beta(α, β)
     */
    inline double GetLogBetaFunction() const { return B.GetLogBetaFunction(); }
};

#endif // BETAPRIMERAND_H

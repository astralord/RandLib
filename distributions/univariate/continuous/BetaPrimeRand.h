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
class RANDLIBSHARED_EXPORT BetaPrimeRand : public BetaRand
{
public:
    BetaPrimeRand(double shape1 = 1, double shape2 = 1);
    std::string Name() const override;
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

    std::complex<double> CF(double t) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};

#endif // BETAPRIMERAND_H

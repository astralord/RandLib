#ifndef BETAPRIMERAND_H
#define BETAPRIMERAND_H

#include "BetaRand.h"

/**
 * @brief The BetaPrimeRand class
 *
 * f(x|\alpha, \beta) = x^{\alpha-1} (1+x)^{-\alpha - \beta} / B(\alpha, \beta)
 *
 * X ~ Beta-Prime(\alpha, \beta)
 * If Y ~ Beta(\alpha, \beta), then X = Y / (1 - Y)
 * If Y ~ Gamma(\alpha), Z ~ Gamma(\beta), then X = Y / (Z - Y)
 */
class RANDLIBSHARED_EXPORT BetaPrimeRand : public BetaRand
{
public:
    BetaPrimeRand(double shape1 = 1, double shape2 = 1);
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return SEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    void sample(std::vector<double> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // BETAPRIMERAND_H

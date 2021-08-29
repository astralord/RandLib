#ifndef BETABINOMIALRAND_H
#define BETABINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/BetaRand.h"
#include <functional>

/**
 * @brief The BetaBinomialRand class <BR>
 * Beta-Binomial distribution
 *
 * Notation: X ~ BB(n, α, β)
 *
 * Related distributions: <BR>
 * If X ~ Binomial(n, p), where p ~ Beta(α, β), then X ~ BB(n, α, β)
 */
template< typename IntType = int >
class RANDLIBSHARED_EXPORT BetaBinomialRand : public DiscreteDistribution<IntType>
{
    IntType n = 1; ///< number of experiments
    double pmfCoef = 0; ///< log(n!) - log(Γ(α + β + n)) - log(B(α, β))
    BetaRand<double> B{};

public:
    BetaBinomialRand(IntType number = 1, double shape1 = 1, double shape2 = 1);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    IntType MinValue() const override { return 0; }
    IntType MaxValue() const override { return n; }

    void SetParameters(IntType number, double shape1, double shape2);
    inline IntType GetNumber() const { return n; }
    inline double GetAlpha() const { return B.GetAlpha(); }
    inline double GetBeta() const { return B.GetBeta(); }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;

private:
    IntType VariateUniform() const;
    IntType VariateBeta() const;

public:
    IntType Variate() const override;
    void Sample(std::vector<IntType> &outputData) const override;

    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;
};

#endif // BETABINOMIALRAND_H

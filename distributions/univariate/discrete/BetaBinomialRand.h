#ifndef BETABINOMIALRAND_H
#define BETABINOMIALRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/BetaRand.h"
#include <functional>

/**
 * @brief The BetaBinomialRand class
 */
class RANDLIBSHARED_EXPORT BetaBinomialRand : public DiscreteDistribution
{
    int n;
    BetaRand B;

public:
    BetaBinomialRand(int number, double shape1, double shape2);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return n; }

    void SetParameters(int number, double shape1, double shape2);
    inline int GetNumber() const { return n; }
    inline double GetAlpha() const { return B.GetAlpha(); }
    inline double GetBeta() const { return B.GetBeta(); }

    double P(int k) const override;
    double F(int k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // BETABINOMIALRAND_H

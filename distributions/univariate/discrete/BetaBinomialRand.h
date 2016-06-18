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
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return n; }

    void setParameters(int number, double shape1, double shape2);
    inline int getNumber() const { return n; }
    inline double getAlpha() const { return B.getAlpha(); }
    inline double getBeta() const { return B.getBeta(); }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;

    double Mean() const override;
    double Variance() const override;

    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // BETABINOMIALRAND_H

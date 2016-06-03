#ifndef ZETARAND_H
#define ZETARAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The ZetaRand class
 * Zeta distribution
 * X ~ Zeta(s)
 *
 * P(X = k) = 1 / (k^s * zeta(s))
 */
class RANDLIBSHARED_EXPORT ZetaRand : public DiscreteDistribution
{
    double s, sm1;
    double zetaSInv; /// 1.0 / zeta(s)
    double b;
public:
    explicit ZetaRand(double exponent = 1.0);
    std::string name() override;

    void setExponent(double exponent);
    inline double getExponent() { return s; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;

    inline double getInverseZetaFunction() { return zetaSInv; }
};

#endif // ZETARAND_H

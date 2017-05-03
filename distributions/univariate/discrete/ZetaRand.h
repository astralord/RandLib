#ifndef ZETARAND_H
#define ZETARAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The ZetaRand class
 * Zeta distribution
 *
 * P(X = k) = 1 / (k^s * ζ(s))
 *
 * Notation: X ~ Zeta(s)
 */
class RANDLIBSHARED_EXPORT ZetaRand : public DiscreteDistribution
{
    double s, sm1;
    /// 1.0 / ζ(s)
    double zetaSInv;
    double b;
public:
    explicit ZetaRand(double exponent = 1.0);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 1; }
    int MaxValue() const override { return INT_MAX; }

    void SetExponent(double exponent);
    inline double GetExponent() const { return s; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double GetInverseZetaFunction() const { return zetaSInv; }
};

#endif // ZETARAND_H

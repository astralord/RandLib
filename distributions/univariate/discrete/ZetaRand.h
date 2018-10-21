#ifndef ZETARAND_H
#define ZETARAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The ZetaRand class <BR>
 * Zeta distribution
 *
 * P(X = k) = 1 / (k^s * ζ(s))
 *
 * Notation: X ~ Zeta(s)
 */
template < typename IntType = int >
class RANDLIBSHARED_EXPORT ZetaRand : public DiscreteDistribution<IntType>
{
    double s = 2; ///< exponent
    double sm1 = 1; ///< s - 1
    double zetaS = M_PI_SQ / 6.0; ///< ζ(s), where ζ stands for Riemann zeta-function
    double logZetaS = 2 * M_LNPI - M_LN2 - M_LN3;///< ln(ζ(s))
    double b = 0.5; ///< 1 - 2^(1-s)

public:
    explicit ZetaRand(double exponent = 2.0);
    String Name() const override;

    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    IntType MinValue() const override { return 1; }
    IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetExponent(double exponent);
    inline double GetExponent() const { return s; }

    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    IntType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;
    long double Moment(int n) const;
    long double ThirdMoment() const override { return Moment(3); }
    long double FourthMoment() const override { return Moment(4); }

    inline long double GetZetaFunction() const { return zetaS; }
    inline long double GetLogZetaFunction() const { return logZetaS; }
};

#endif // ZETARAND_H

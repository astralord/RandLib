#ifndef LOGARITHMICRAND_H
#define LOGARITHMICRAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The LogarithmicRand class <BR>
 * Logarithmic distribution
 *
 * P(X = k) = -p^k / [k log(1 - p)]
 *
 * X ~ Log(p)
 */
template < typename IntType = int >
class RANDLIBSHARED_EXPORT LogarithmicRand : public DiscreteDistribution<IntType>
{
    double p = 0.5; ///< parameter of distribution
    double logProb = -M_LN2; ///< log(p)
    double log1mProb = -M_LN2; ///< log(q)
public:
    explicit LogarithmicRand(double probability);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    IntType MinValue() const override { return 1; }
    IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetProbability(double probability);
    inline double GetProbability() const { return p; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
private:
    /**
     * @fn betaFun
     * @param a
     * @return B(p, a, 0), where B(x, a, b) denotes incomplete beta function,
     * using series expansion (converges for x < 1)
     */
    double betaFun(IntType a) const;
public:
    double F(const IntType & k) const override;
    double S(const IntType & k) const override;
    IntType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // LOGARITHMICRAND_H

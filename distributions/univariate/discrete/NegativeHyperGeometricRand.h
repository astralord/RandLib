#ifndef NEGATIVEHYPERGEOMETRICRAND_H
#define NEGATIVEHYPERGEOMETRICRAND_H

#include "DiscreteDistribution.h"
#include "BernoulliRand.h"

/**
 * @brief The NegativeHyperGeometricRand class <BR>
 * Negative hypergeometric distribution
 *
 * Notation: X ~ NHG(N, M, m)
 */
template < typename IntType = int >
class RANDLIBSHARED_EXPORT NegativeHyperGeometricRand : public DiscreteDistribution<IntType>
{
    IntType N = 1; ///< size of population
    IntType M = 1; ///< total amount of successes
    IntType m = 1; ///< limiting number of successes
    double pmfCoef = 1; ///< C(N, M)
    double p0 = 1; ///< M / N

public:
    NegativeHyperGeometricRand(IntType totalSize = 1, IntType totalSuccessesNum = 1, IntType limitSuccessesNum = 1);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    IntType MinValue() const override { return 0; }
    IntType MaxValue() const override { return N - M; }

    void SetParameters(IntType totalSize, IntType totalSuccessesNum, IntType limitSuccessesNum);
    inline IntType GetTotalSize() { return N; }
    inline IntType GetTotalSuccessesNum() { return M; }
    inline IntType GetLimitSuccessesNum() { return m; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    IntType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
};

#endif // NEGATIVEHYPERGEOMETRICRAND_H

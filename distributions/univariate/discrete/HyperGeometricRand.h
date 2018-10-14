#ifndef HYPERGEOMETRICRAND_H
#define HYPERGEOMETRICRAND_H

#include "DiscreteDistribution.h"
#include "BernoulliRand.h"
#include "BetaBinomialRand.h"

/**
 * @brief The HyperGeometricRand class <BR>
 * Hypergeometric distribution
 *
 * X ~ HG(N, K, n)
 */
template < typename IntType = int >
class RANDLIBSHARED_EXPORT HyperGeometricRand : public DiscreteDistribution<IntType>
{
    IntType N = 1; ///< population size
    IntType K = 1; /// number of possible successes
    IntType n = 1; /// number of draws
    double pmfCoef = 1; ///< C(N, n)
    double p0 = 1; ///< K/N

public:
    HyperGeometricRand(IntType totalSize = 1, IntType drawsNum = 1, IntType successesNum = 1);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    IntType MinValue() const override { return std::max(static_cast<IntType>(0), n - N + K); }
    IntType MaxValue() const override { return std::min(n, K); }

    void SetParameters(IntType totalSize, IntType drawsNum, IntType successesNum);
    inline IntType GetTotalSize() { return N; }
    inline IntType GetDrawsNum() { return n; }
    inline IntType GetSuccessesNum() { return K; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    IntType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;
};

#endif // HYPERGEOMETRICRAND_H

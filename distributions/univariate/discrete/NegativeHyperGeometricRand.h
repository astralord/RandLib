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
class RANDLIBSHARED_EXPORT NegativeHyperGeometricRand : public DiscreteDistribution
{
    int N = 1; ///< size of population
    int M = 1; ///< total amount of successes
    int m = 1; ///< limiting number of successes
    double pmfCoef = 1; ///< C(N, M)
    double p0 = 1; ///< M / N

public:
    NegativeHyperGeometricRand(size_t totalSize, size_t totalSuccessesNum, size_t limitSuccessesNum);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return N - M; }

    void SetParameters(size_t totalSize, size_t totalSuccessesNum, size_t limitSuccessesNum);
    inline size_t GetTotalSize() { return N; }
    inline size_t GetTotalSuccessesNum() { return M; }
    inline size_t GetLimitSuccessesNum() { return m; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
};

#endif // NEGATIVEHYPERGEOMETRICRAND_H

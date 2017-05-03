#ifndef NEGATIVEHYPERGEOMETRICRAND_H
#define NEGATIVEHYPERGEOMETRICRAND_H

#include "DiscreteDistribution.h"
#include "BernoulliRand.h"

/**
 * @brief The NegativeHyperGeometricRand class
 * Negative hypergeometric distribution
 *
 * Notation: X ~ NHG(N, M, m)
 */
class RANDLIBSHARED_EXPORT NegativeHyperGeometricRand : public DiscreteDistribution
{
    int N, M, m;
    /// C(N, M)
    double pmfCoef;
    double p0;

public:
    NegativeHyperGeometricRand(int totalSize, int totalSuccessesNum, int limitSuccessesNum);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return N - M; }

    void SetParameters(int totalSize, int totalSuccessesNum, int limitSuccessesNum);
    inline int GetTotalSize() { return N; }
    inline int GetTotalSuccessesNum() { return M; }
    inline int GetLimitSuccessesNum() { return m; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
};

#endif // NEGATIVEHYPERGEOMETRICRAND_H

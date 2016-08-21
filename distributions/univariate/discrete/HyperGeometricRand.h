#ifndef HYPERGEOMETRICRAND_H
#define HYPERGEOMETRICRAND_H

#include "DiscreteDistribution.h"
#include "BernoulliRand.h"

/**
 * @brief The HyperGeometricRand class
 * Hyper-geometric distribution
 */
class RANDLIBSHARED_EXPORT HyperGeometricRand : public DiscreteDistribution
{
    int N, K, n;
    double pdfDenominator; /// C(N, n)
    double p0;

public:
    HyperGeometricRand(int totalSize, int drawsNum, int successesNum);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return std::max(0, n - N + K); }
    int MaxValue() const override { return std::min(n, K); }

    void SetParameters(int totalSize, int drawsNum, int successesNum);
    inline int GetTotalSize() { return N; }
    inline int GetDrawsNum() { return n; }
    inline int GetSuccessesNum() { return K; }

    double P(int k) const override;
    double F(int k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // HYPERGEOMETRICRAND_H

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
    std::string name() override;

    void setParameters(int totalSize, int drawsNum, int successesNum);
    inline int getTotalSize() { return N; }
    inline int getDrawsNum() { return n; }
    inline int getSuccessesNum() { return K; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // HYPERGEOMETRICRAND_H

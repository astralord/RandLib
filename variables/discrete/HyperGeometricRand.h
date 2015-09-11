#ifndef HYPERGEOMETRICRAND_H
#define HYPERGEOMETRICRAND_H

#include "DiscreteRand.h"
#include "BernoulliRand.h"

/**
 * @brief The HyperGeometricRand class
 */
class RANDLIBSHARED_EXPORT HyperGeometricRand : public DiscreteRand<int>
{
    size_t N, K, n;
    double pdfDenominator; /// C(N, n)
public:
    HyperGeometricRand(size_t totalSize, size_t drawsNum, size_t successesNum);
    virtual std::string name() override;

    void setParameters(size_t totalSize, size_t drawsNum, size_t successesNum);
    inline size_t getTotalSize() { return N; }
    inline size_t getDrawsNum() { return n; }
    inline size_t getSuccessesNum() { return K; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override;
    double Var() const override;

    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // HYPERGEOMETRICRAND_H

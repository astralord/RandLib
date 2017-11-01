#ifndef UNIFORM_DISCRETE_RAND_H
#define UNIFORM_DISCRETE_RAND_H

#include "DiscreteDistribution.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The UniformDiscreteRand class <BR>
 * Uniform discrete distribution
 *
 * Notation: X ~ U(a, b)
 *
 * P(X = k) = 1 / (b - a + 1) for a <= k <= b
 */
class RANDLIBSHARED_EXPORT UniformDiscreteRand : public DiscreteDistribution
{
    size_t n = 1; ///< number of possible outcomes
    int a = 0; ///< min bound
    int b = 0; ///< max bound
    double nInv = 1; ///< 1/n
    double logN = 0; ///< log(n)
    unsigned long long MAX_RAND_UNBIASED = RandGenerator::MaxValue();///< constant for unbiased generator

public:
    UniformDiscreteRand(int minValue = 0, int maxValue = 1);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return a; }
    int MaxValue() const override { return b; }

    void SetBoundaries(int minValue, int maxValue);

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;
    double LikelihoodFunction(const std::vector<int> &sample) const override;
    double LogLikelihoodFunction(const std::vector<int> &sample) const override;
};

#endif // UNIFORM_DISCRETE_RAND_H

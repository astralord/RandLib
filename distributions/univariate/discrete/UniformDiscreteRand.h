#ifndef UNIFORM_DISCRETE_RAND_H
#define UNIFORM_DISCRETE_RAND_H

#include "DiscreteDistribution.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The UniformDiscreteRand class
 * Uniform discrete distribution
 *
 *  X ~ U(a, b)
 *
 * P(X = k) = 1_{k \in [a, b]} / (b - a + 1)
 */
class RANDLIBSHARED_EXPORT UniformDiscreteRand : public DiscreteDistribution
{
    int n, a, b;
    double nInv; /// 1.0 / n
    double logN; /// log(n)

public:
    UniformDiscreteRand(int minValue = 0, int maxValue = 1);
    std::string Name() const override;
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
    double Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;
    double Likelihood(const std::vector<int> &sample) const override;
    double LogLikelihood(const std::vector<int> &sample) const override;
};

#endif // UNIFORM_DISCRETE_RAND_H

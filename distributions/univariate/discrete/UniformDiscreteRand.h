#ifndef UNIFORM_DISCRETE_RAND_H
#define UNIFORM_DISCRETE_RAND_H

#include "DiscreteDistribution.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The UniformDiscreteRand class
 * Uniform discrete distribution
 * X ~ U(a, b)
 *
 * P(X = k) = 1_{x \in [a, b]} / (b - a + 1)
 */
class RANDLIBSHARED_EXPORT UniformDiscreteRand : public DiscreteDistribution
{
    int n, a, b;
    double nInv; /// 1.0 / n

public:
    UniformDiscreteRand(int minValue = 0, int maxValue = 1);
    std::string name() override;

    void setBoundaries(int minValue, int maxValue);
    inline int getMinValue() const { return a; }
    inline int getMaxValue() const { return b; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;

    inline double Entropy() { return std::log(n); }
};

#endif // UNIFORM_DISCRETE_RAND_H

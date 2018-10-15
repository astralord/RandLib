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
template< typename IntType = int >
class RANDLIBSHARED_EXPORT UniformDiscreteRand : public DiscreteDistribution<IntType>
{
    size_t n = 1; ///< number of possible outcomes
    IntType a = 0; ///< min bound
    IntType b = 0; ///< max bound
    double nInv = 1; ///< 1/n
    double logN = 0; ///< log(n)
    unsigned long long MAX_RAND_UNBIASED = this->localRandGenerator.MaxValue();///< constant for unbiased generator

public:
    UniformDiscreteRand(IntType minValue = 0, IntType maxValue = 1);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    IntType MinValue() const override { return a; }
    IntType MaxValue() const override { return b; }

    void SetBoundaries(IntType minValue, IntType maxValue);

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;

    IntType Variate() const override;
    static IntType StandardVariate(IntType minValue = 0, IntType maxValue = 1, RandGenerator &randGenerator = ProbabilityDistribution<IntType>::staticRandGenerator);

    long double Mean() const override;
    long double Variance() const override;
    IntType Median() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    long double Entropy() const override;
    double LikelihoodFunction(const std::vector<IntType> &sample) const override;
    double LogLikelihoodFunction(const std::vector<IntType> &sample) const override;

    /**
     * @fn Fit
     * fit bounds via maximum-likelihood method
     * @param sample
     */
    void Fit(const std::vector<IntType> &sample);
};

#endif // UNIFORM_DISCRETE_RAND_H

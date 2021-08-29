#ifndef CATEGORICALRAND_H
#define CATEGORICALRAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The CategoricalRand class <BR>
 *
 * P(X = k) = p_k for k = {0, ..., K-1}
 *
 * Notation: X ~ Cat(p_0, ..., p_{K-1})
 *
 * Related distributions: <BR>
 * X ~ Multin(1, p_0, ..., p_{K-1}) <BR>
 * If X ~ Bernoulli(p), then X ~ Cat(1 - p, p) <BR>
 * If X ~ Uniform-Discrete(0, K), then X ~ Cat(p, ..., p) with p = 1 / (K + 1)
 */
template < typename IntType = int >
class RANDLIBSHARED_EXPORT CategoricalRand : public DiscreteDistribution<IntType>
{
    std::vector<double> prob{1.0}; ///< vector of probabilities
    IntType K = 1; ///< number of possible outcomes

public:
    explicit CategoricalRand(std::vector<double>&& probabilities = {1.0});
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    IntType MinValue() const override { return 0; }
    IntType MaxValue() const override { return K - 1; }

    void SetProbabilities(std::vector<double>&& probabilities);
    std::vector<double> GetProbabilities() { return prob; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    IntType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;

private:
    IntType quantileImpl(double p) const override;
    IntType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};


#endif // CATEGORICALRAND_H

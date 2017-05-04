#ifndef CATEGORICALRAND_H
#define CATEGORICALRAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The CategoricalRand class <BR>
 *
 * P(X = k) = p_k for k = {0, ..., K}
 *
 * Notation: X ~ Cat(p_0, p_1, ..., p_K, K)
 *
 * Related distributions:
 * X ~ Multin(1, p_0, p_1, ..., p_K, K)
 * If X ~ Bernoulli(p), then X ~ Cat(1-p, p, 2)
 * If X ~ Uniform-Discrete(0, K), then X ~ Cat(p, ..., p, K) with p = 1 / (K + 1)
 */
class RANDLIBSHARED_EXPORT CategoricalRand : public DiscreteDistribution
{
    std::vector<double> prob;
    int K;
    double q; ///< P(X = 0)

public:
    explicit CategoricalRand(std::vector<double>&& probabilities);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 0; }
    int MaxValue() const override { return K; }

    void SetProbabilities(std::vector<double>&& probabilities);
    std::vector<double> GetProbabilities() { return prob; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;

private:
    int quantileImpl(double p) const override;
    int quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};


#endif // CATEGORICALRAND_H

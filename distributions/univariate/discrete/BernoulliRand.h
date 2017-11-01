#ifndef BERNOULLIRAND_H
#define BERNOULLIRAND_H

#include "BinomialRand.h"

/**
 * @brief The BernoulliRand class <BR>
 * Bernoulli distribution
 *
 * P(X = k) = p * 1_{k = 1} + (1 - p) * 1_{k = 0}
 *
 * Notation: X ~ Bernoulli(p)
 *
 * Related distributions: <BR>
 * X ~ Binomial(1, p) <BR>
 * X ~ Multin(1, 1 - p, p) <BR>
 * 2X - 1 ~ Rademacher
 */
class RANDLIBSHARED_EXPORT BernoulliRand : public BinomialDistribution
{
    unsigned long long boundary = 0;///< coefficient for faster random number generation

public:
    explicit BernoulliRand(double probability = 0.5);
    String Name() const override;

public:
    void SetProbability(double probability);

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;
    static int Variate(double probability);
    static int StandardVariate();
    void Sample(std::vector<int> &outputData) const override;

    inline double Entropy();
};

#endif // BERNOULLIRAND_H

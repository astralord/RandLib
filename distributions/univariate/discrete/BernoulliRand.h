#ifndef BERNOULLIRAND_H
#define BERNOULLIRAND_H

#include "BinomialRand.h"

/**
 * @brief The BernoulliRand class
 * Bernoulli distribution
 *
 * P(X = k) = p * 1_{k = 1} + (1 - p) * 1_{k = 0}
 *
 * Notation: X ~ Bernoulli(p)
 *
 * Related distributions:
 * X ~ Binomial(1, p)
 * 2X - 1 ~ Rademacher
 */
class RANDLIBSHARED_EXPORT BernoulliRand : public BinomialRand
{
    unsigned long long boundary;

public:
    explicit BernoulliRand(double probability = 0.5);
    std::string Name() const override;

protected:
    using BinomialRand::SetParameters;

public:
    void SetProbability(double probability);

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;
    static int Variate(double p);
    static int StandardVariate();

    inline double Entropy();
};

#endif // BERNOULLIRAND_H

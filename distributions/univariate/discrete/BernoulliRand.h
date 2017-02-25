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

    double P(int k) const override;
    double logP(int k) const override;
    double F(int k) const override;
    double S(int k) const override;
    int Variate() const override;
    static int Variate(double p);
    static int StandardVariate();

    inline double Entropy();

    bool FitMLE(const std::vector<int> &sample);
    bool FitMM(const std::vector<int> &sample);
    bool FitBayes(const std::vector<int> &sample, BetaRand &priorDistribution);
};

#endif // BERNOULLIRAND_H

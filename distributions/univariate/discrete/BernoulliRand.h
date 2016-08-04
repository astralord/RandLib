#ifndef BERNOULLIRAND_H
#define BERNOULLIRAND_H

#include "BinomialRand.h"

/**
 * @brief The BernoulliRand class
 * Bernoulli distribution
 * X ~ Bernoulli(p)
 *
 * P(X = 0) = 1 - p
 * P(X = 1) = p
 *
 * X ~ Binomial(1, p)
 */
class RANDLIBSHARED_EXPORT BernoulliRand : public BinomialRand
{
    unsigned long long boundary;

public:
    explicit BernoulliRand(double probability = 0.5);
    std::string name() const override;

protected:
    using BinomialRand::setParameters;

public:
    void setProbability(double probability);

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;
    static int variate(double p);
    static int standardVariate();

    inline double Entropy();

    bool fitMLE(const std::vector<int> &sample);
    bool fitMM(const std::vector<int> &sample);
    bool fitBayes(const std::vector<int> &sample, BetaRand &priorDistribution);
};

#endif // BERNOULLIRAND_H

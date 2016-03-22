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
    std::string name() override;

protected:
    using BinomialRand::setParameters;

public:
    void setProbability(double probability);

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double p);
    static double standardVariate();

    inline double Entropy();
};

#endif // BERNOULLIRAND_H

#ifndef GEOMETRICRAND_H
#define GEOMETRICRAND_H

#include "NegativeBinomialRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The GeometricRand class
 * Geometric distribution
 *
 * Notation: X ~ Geometric(p)
 *
 * P(X = k) = p (1 - p)^k
 *
 * Related distributions:
 * X ~ NB(1, p)
 */
class RANDLIBSHARED_EXPORT GeometricRand : public PascalRand
{
public:
    explicit GeometricRand(double probability = 0.5);
    std::string name() const override;

protected:
    using NegativeBinomialRand::setParameters;

public:
    void setProbability(double probability);

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;
    static int variate(double probability);

    void sample(std::vector<int> &outputData) const override;

    double Median() const override;

    double Entropy() const;

    bool fitMLE(const std::vector<int> &sample);
    bool fitMM(const std::vector<int> &sample);
    bool fitBayes(const std::vector<int> &sample, BetaRand &priorDistribution);
};

#endif // GEOMETRICRAND_H

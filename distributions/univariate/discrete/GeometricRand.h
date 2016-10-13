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
    std::string Name() const override;

protected:
    using NegativeBinomialRand::SetParameters;

public:
    void SetProbability(double probability);

    double P(int k) const override;
    double F(int k) const override;
    int Variate() const override;
    static int Variate(double probability);

    void Sample(std::vector<int> &outputData) const override;

    double Median() const override;

    double Entropy() const;

    bool FitMLE(const std::vector<int> &sample);
    bool FitMM(const std::vector<int> &sample);
    bool FitBayes(const std::vector<int> &sample, BetaRand &priorDistribution);
};

#endif // GEOMETRICRAND_H

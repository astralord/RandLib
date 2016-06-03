#ifndef GEOMETRICRAND_H
#define GEOMETRICRAND_H

#include "NegativeBinomialRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The GeometricRand class
 * Geometric distribution
 * X ~ Geometric(p)
 *
 * P(X = k) = p (1 - p)^k
 *
 * X ~ NB(1, p)
 */
class RANDLIBSHARED_EXPORT GeometricRand : public PascalRand
{
public:
    explicit GeometricRand(double probability = 0.5);
    std::string name() override;

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
    
    bool checkValidity(const std::vector<double> &sample);

    bool fitProbabilityMLE(const std::vector<double> &sample);
    bool fitProbabilityMM(const std::vector<double> &sample);
    bool fitProbabilityBayes(const std::vector<double> &sample, BetaRand &priorDistribution);
};

#endif // GEOMETRICRAND_H

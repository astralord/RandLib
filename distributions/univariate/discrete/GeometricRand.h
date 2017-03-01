#ifndef GEOMETRICRAND_H
#define GEOMETRICRAND_H

#include "NegativeBinomialRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The GeometricRand class
 * Geometric distribution
 *
 * P(X = k) = p (1 - p)^k
 *
 * Notation: X ~ Geometric(p)
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
    using NegativeBinomialRand::FitMM;
    using NegativeBinomialRand::FitNumberMM;

public:
    void SetProbability(double probability);

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;
    static int Variate(double probability);

    void Sample(std::vector<int> &outputData) const override;

    double Median() const override;

    double Entropy() const;

    void FitMLE(const std::vector<int> &sample);
    void FitMM(const std::vector<int> &sample);
    BetaRand FitBayes(const std::vector<int> &sample, const BetaRand &priorDistribution);
};

#endif // GEOMETRICRAND_H

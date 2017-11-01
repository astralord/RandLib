#ifndef GEOMETRICRAND_H
#define GEOMETRICRAND_H

#include "NegativeBinomialRand.h"

/**
 * @brief The GeometricRand class <BR>
 * Geometric distribution
 *
 * P(X = k) = p (1 - p)^k
 *
 * Notation: X ~ Geometric(p)
 *
 * Related distributions: <BR>
 * X ~ NB(1, p)
 */
class RANDLIBSHARED_EXPORT GeometricRand : public NegativeBinomialDistribution<int>
{
public:
    explicit GeometricRand(double probability = 0.5) : NegativeBinomialDistribution<int>(1, probability) {}
    String Name() const override;

public:
    void SetProbability(double probability);

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;
    static int Variate(double probability);

    void Sample(std::vector<int> &outputData) const override;

    int Median() const override;

    double Entropy() const;
};

#endif // GEOMETRICRAND_H

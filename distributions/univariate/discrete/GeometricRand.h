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
template < typename IntType = int >
class RANDLIBSHARED_EXPORT GeometricRand : public PascalRand<IntType>
{
public:
    explicit GeometricRand(double probability = 0.5) : PascalRand<IntType>(1, probability) {}
    String Name() const override;

public:
    void SetProbability(double probability);

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    double S(const IntType & k) const override;
    IntType Variate() const override;
    static IntType Variate(double probability, RandGenerator &randGenerator = ProbabilityDistribution<IntType>::staticRandGenerator);

    void Sample(std::vector<IntType> &outputData) const override;

    IntType Median() const override;

    double Entropy() const;
};

#endif // GEOMETRICRAND_H

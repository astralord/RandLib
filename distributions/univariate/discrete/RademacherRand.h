#ifndef RADEMACHERRAND_H
#define RADEMACHERRAND_H

#include "DiscreteDistribution.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The RademacherRand class <BR>
 * Rademacher distribution
 *
 * P(X = k) = 0.5 * 1_{|k| = 1}
 *
 * Notation: X ~ Rademacher
 *
 * Related distributions: <BR>
 * If Y ~ Bernoulli(0.5), then 2Y - 1 ~ Rademacher
 */
class RANDLIBSHARED_EXPORT RademacherRand : public DiscreteDistribution
{
public:
    RademacherRand();
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return -1; }
    int MaxValue() const override { return 1; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    int quantileImpl(double p) const override;
    int quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy();
};

#endif // RADEMACHERRAND_H

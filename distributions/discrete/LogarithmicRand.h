#ifndef LOGARITHMICRAND_H
#define LOGARITHMICRAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The LogarithmicRand class
 * Logarithmic distribution
 * X ~ Logarithmic(p)
 *
 * P(X = k) = -p^k / [k \log(1 - p)];
 */
class RANDLIBSHARED_EXPORT LogarithmicRand : public DiscreteDistribution
{
    double p, q;
    double logQInv; /// 1 / log(q)
public:
    explicit LogarithmicRand(double probability);
    std::string name() override;

    void setProbability(double probability);
    inline double getProbability() const { return p; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Mode() const override;
};

#endif // LOGARITHMICRAND_H

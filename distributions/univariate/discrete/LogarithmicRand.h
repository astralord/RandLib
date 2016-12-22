#ifndef LOGARITHMICRAND_H
#define LOGARITHMICRAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The LogarithmicRand class
 * Logarithmic distribution
 *
 * P(X = k) = -p^k / [k log(1 - p)]
 *
 * X ~ Logarithmic(p)
 */
class RANDLIBSHARED_EXPORT LogarithmicRand : public DiscreteDistribution
{
    double p, q;
    double logQInv; /// 1 / log(q)
public:
    explicit LogarithmicRand(double probability);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 1; }
    int MaxValue() const override { return INT_MAX; }

    void SetProbability(double probability);
    inline double GetProbability() const { return p; }

    double P(int k) const override;
    double F(int k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;

    std::complex<double> CF(double t) const override;
};

#endif // LOGARITHMICRAND_H

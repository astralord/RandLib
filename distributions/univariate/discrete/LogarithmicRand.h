#ifndef LOGARITHMICRAND_H
#define LOGARITHMICRAND_H

#include "DiscreteDistribution.h"

/**
 * @brief The LogarithmicRand class <BR>
 * Logarithmic distribution
 *
 * P(X = k) = -p^k / [k log(1 - p)]
 *
 * X ~ Log(p)
 */
class RANDLIBSHARED_EXPORT LogarithmicRand : public DiscreteDistribution
{
    double p = 1; ///< parameter of distribution
    double logProb = 0; ///< log(p)
    double log1mProb = -INFINITY; ///< log(q)
public:
    explicit LogarithmicRand(double probability);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 1; }
    int MaxValue() const override { return INT_MAX; }

    void SetProbability(double probability);
    inline double GetProbability() const { return p; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
private:
    /**
     * @brief betaFun
     * @param a
     * @return B(p, a, 0), where B(x, a, b) denotes incomplete beta function,
     * using series expansion (converges for x < 1)
     */
    double betaFun(int a) const;
public:
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // LOGARITHMICRAND_H

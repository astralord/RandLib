#ifndef TRINOMIALRAND_H
#define TRINOMIALRAND_H

#include "DiscreteBivariateDistribution.h"
#include "../univariate/discrete/BinomialRand.h"

/**
 * @brief The TrinomialRand class <BR>
 * Trinomial distribution
 *
 * Notation: X ~ Trin(n, p_1, p_2)
 *
 * Related distributions: <BR>
 * X ~ Multin(n, 1 - p_1 - p_2, p_1, p_2)
 */
class RANDLIBSHARED_EXPORT TrinomialRand : public DiscreteBivariateDistribution<BinomialRand, BinomialRand>
{
    int n = 1; ///< number of trials
    double log1mProb = -M_LN3; ///< log(1 - p_1 - p_2)
    double p1_1mp2 = 0.5; ///< p_1 / (1 - p_2)
    double p2_1mp1 = 0.5; ///< p_2 / (1 - p_1)
public:
    TrinomialRand(int number, double probability1, double probability2);
    String Name() const override;

    void SetParameters(int number, double probability1, double probability2);
    inline int GetNumber() const { return n; }
    inline double GetFirstProbability() const { return X.GetProbability(); }
    inline double GetSecondProbability() const { return Y.GetProbability(); }

    double P(const IntPair &point) const override;
    double logP(const IntPair &point) const override;
    double F(const IntPair & point) const override;
    IntPair Variate() const override;

    double Correlation() const override;
};

#endif // TRINOMIALRAND_H

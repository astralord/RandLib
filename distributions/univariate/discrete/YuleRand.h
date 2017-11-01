#ifndef YULERAND_H
#define YULERAND_H

#include "DiscreteDistribution.h"
#include "GeometricRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/ParetoRand.h"

/**
 * @brief The YuleRand class <BR>
 * Yule distribution
 *
 * Notation: X ~ Yule(ρ)
 *
 * Related distributions: <BR>
 * If Y ~ Pareto(ρ, 1) and Z ~ Geometric(1 / Y), then Z + 1 ~ Yule(ρ)
 */
class RANDLIBSHARED_EXPORT YuleRand : public DiscreteDistribution
{
    double ro = 0; ///< shape ρ
    double lgamma1pRo = 0; /// log(Γ(1 + ρ))
    
    ParetoRand X;
public:
    explicit YuleRand(double shape);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 1; }
    int MaxValue() const override { return INT_MAX; }

    void SetShape(double shape);
    inline double GetShape() const { return ro; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;
    static int Variate(double shape);

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // YULERAND_H

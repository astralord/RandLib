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
template < typename IntType = int >
class RANDLIBSHARED_EXPORT YuleRand : public DiscreteDistribution<IntType>
{
    double rho = 0; ///< shape ρ
    double lgamma1pRo = 0; /// log(Γ(1 + ρ))
    
    ParetoRand<double> X;
public:
    explicit YuleRand(double shape);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    IntType MinValue() const override { return 1; }
    IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetShape(double shape);
    inline double GetShape() const { return rho; }

    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    double S(const IntType & k) const override;
    IntType Variate() const override;
    static IntType Variate(double shape, RandGenerator &randGenerator = ProbabilityDistribution<IntType>::staticRandGenerator);
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;
};

#endif // YULERAND_H

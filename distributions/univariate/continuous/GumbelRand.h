#ifndef GUMBELRAND_H
#define GUMBELRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GumbelRand class <BR>
 * Gumbel distribution
 *
 * Notation: X ~ Gumbel(μ, β)
 *
 * Related distributions: <BR>
 * exp(-(X - μ) / β) ~ Exp(1)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT GumbelRand : public ContinuousDistribution<RealType>
{
    double mu = 0; ///< location μ
    double beta = 1; ///< scale β
    double logBeta = 0; ///< log(β)
public:
    GumbelRand(double location = 0, double scale = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    RealType MinValue() const override { return -INFINITY; }
    RealType MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return beta; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;
    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

public:
    double Entropy() const;
};

#endif // GUMBELRAND_H

#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "StableRand.h"

/**
 * @brief The CauchyRand class <BR>
 * Cauchy distribution
 *
 * f(x | μ, σ) = σ / [π (σ^2 + (x - μ)^2)]
 *
 * Notation: X ~ Cauchy(μ, σ)
 *
 * Related distributions: <BR>
 * If X ~ Cauchy(0, 1), then μ + σ * X ~ Cauchy(μ, σ) <BR>
 * X ~ S(1, 0, σ, μ) <BR>
 * If X, Y ~ Normal(0, 1), then X / Y ~ Cauchy(0, 1)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT CauchyRand : public StableDistribution<RealType>
{
public:
    CauchyRand(double location = 0, double scale = 1);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    RealType MinValue() const override { return -INFINITY; }
    RealType MaxValue() const override { return INFINITY; }

public:
    double f(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;
};

#endif // CAUCHYRAND_H

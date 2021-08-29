#ifndef LEVYRAND_H
#define LEVYRAND_H

#include "StableRand.h"

/**
 * @brief The LevyRand class <BR>
 * Levy distribution
 *
 * f(x | μ, σ) = ((σ exp(σ / (μ - x)) / (2 π (x - μ)^3))^(1/2)
 *
 * Notation: X ~ Levy(μ, σ)
 *
 * Related distributions: <BR>
 * If X ~ Levy(0, 1), then μ + σ * X ~ Levy(μ, σ) <BR>
 * X ~ S(0.5, 1, σ, μ) <BR>
 * If Y ~ Normal(0, 1), then 1 / X^2 ~ Levy(0, 1) <BR>
 * If X ~ Levy(0, σ), then X ~ Inv-Γ(1/2, σ/2)
 */
template < typename RealType = double>
class RANDLIBSHARED_EXPORT LevyRand : public StableDistribution<RealType>
{
public:
    LevyRand(double location = 0, double scale = 1);
    String Name() const override;

public:
    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @fn FitScale
     * Fit scale using maximum-likelihoood estimator
     * @param sample
     */
    void FitScale(const std::vector<RealType> &sample);
};

#endif // LEVYRAND_H

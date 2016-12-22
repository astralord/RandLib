#ifndef LEVYRAND_H
#define LEVYRAND_H

#include "StableRand.h"

/**
 * @brief The LevyRand class
 * Levy distribution
 *
 * f(x | μ, σ) = ((σ exp(σ / (μ - x)) / (2 π (x - μ)^3))^(1/2)
 *
 * Notation: X ~ Levy(μ, σ)
 *
 * Related distributions:
 * X ~ Stable(0.5, 1, σ, μ)
 * If Y ~ Normal(0, 1), then μ + σ / X^2 ~ Levy(μ, σ)
 */
class RANDLIBSHARED_EXPORT LevyRand : public StableRand
{
public:
    LevyRand(double location = 0, double scale = 1);
    std::string Name() const override;

private:
    using StableRand::SetParameters;
    using StableRand::GetExponent;
    using StableRand::GetSkewness;

public:
    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

    static double Variate(double location, double scale);
    static double StandardVariate();
    
    std::complex<double> CF(double t) const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    /// Maximum likelihood estimators
    bool FitScaleMLE(const std::vector<double> &sample);
};

#endif // LEVYRAND_H

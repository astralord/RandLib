#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "StableRand.h"

/**
 * @brief The CauchyRand class
 * Cauchy distribution
 *
 * f(x | μ, σ) = σ / [π (σ^2 + (x - μ)^2)]
 *
 * Notation: X ~ Cauchy(μ, σ)
 *
 * Related distributions:
 * If X ~ Cauchy(0, 1), then μ + σ * X ~ Cauchy(μ, σ)
 * X ~ Stable(1, 0, σ, μ)
 * If X, Y ~ Normal(0, 1), then X/Y ~ Cauchy(0, 1)
 */
class RANDLIBSHARED_EXPORT CauchyRand : public StableRand
{
public:
    CauchyRand(double location = 0, double scale = 1);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

private:
    using StableRand::SetParameters;

public:
    double f(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;

    static double Variate(double location, double scale);
    static double StandardVariate();

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;
};

#endif // CAUCHYRAND_H

#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "StableRand.h"

/**
 * @brief The CauchyRand class
 * Cauchy distribution
 *
 * f(x|μ, σ) = σ / [π (σ^2 + (x - μ)^2)]
 *
 * Notation: X ~ Cauchy(μ, σ)
 *
 * Related distributions:
 * X ~ Stable(1, 0, σ, μ)
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
    double Entropy() const;

    /// Method of moments
    bool FitMM(const std::vector<double> &) { return false; }
};

#endif // CAUCHYRAND_H

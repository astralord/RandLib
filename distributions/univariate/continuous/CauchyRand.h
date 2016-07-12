#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "StableRand.h"

/**
 * @brief The CauchyRand class
 * Cauchy distribution
 * X ~ Cauchy(μ, σ)
 *
 * f(x|μ, σ) = σ / [π (σ^2 + (x - μ)^2)]
 *
 */
class RANDLIBSHARED_EXPORT CauchyRand : public StableRand
{
public:
    CauchyRand(double location = 0, double scale = 1);
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

private:
    using StableRand::setParameters;
    using StableRand::getExponent;
    using StableRand::getSkewness;

public:
    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double scale);
    static double standardVariate();

    std::complex<double> CF(double t) const override;

    static double standardQuantile(double p);
    static double quantile(double p, double mean, double scale);
    double Quantile(double p) const override;

    double Entropy() const;
};

#endif // CAUCHYRAND_H

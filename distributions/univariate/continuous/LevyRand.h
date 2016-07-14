#ifndef LEVYRAND_H
#define LEVYRAND_H

#include "StableRand.h"

/**
 * @brief The LevyRand class
 * Levy distribution
 *
 * f(x|μ, σ) = ((σ exp(σ / (μ - x)) / (2 π (x - μ)^3))^(1/2)
 *
 * Notation: X ~ Levy(μ, σ)
 *
 * Related distributions:
 * X ~ Stable(0.5, 1, σ, μ)
 */
class RANDLIBSHARED_EXPORT LevyRand : public StableRand
{
public:
    LevyRand(double location = 0, double scale = 1);
    std::string name() const override;

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
    double Quantile(double p) const override;
    
    /// Verify that all elements of sample can have this distribution
    bool checkValidity(const std::vector<double> &sample);
    
    /// Maximum likelihood estimators
    bool fitScaleMLE(const std::vector<double> &sample);
};

#endif // LEVYRAND_H

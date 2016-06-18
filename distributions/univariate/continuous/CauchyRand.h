#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "StableRand.h"

/**
 * @brief The CauchyRand class
 * Cauchy distribution
 * X ~ Cauchy(\mu, \sigma)
 *
 * f(x|\mu, \sigma) = \sigma / [\pi (\sigma^2 + (x - \mu)^2)]
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

    double Quantile(double p) const override;

    double Entropy() const;
};

#endif // CAUCHYRAND_H

#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "StableRand.h"

/**
 * @brief The CauchyRand class
 *
 * f(x|x_0, \gamma) = \gamma / [\pi (\gamma^2 + (x - x_0)^2)]
 *
 * X ~ Cauchy(x_0, \gamma)
 */
class RANDLIBSHARED_EXPORT CauchyRand : public StableRand
{
public:
    CauchyRand(double location = 0, double scale = 1);
    std::string name() override;

private:
    using StableRand::setParameters;

public:
    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double scale);
    static double standardVariate();

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;

    double Entropy() const;
};

#endif // CAUCHYRAND_H

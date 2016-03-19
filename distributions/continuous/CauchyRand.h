#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "ContinuousRand.h"

/**
 * @brief The CauchyRand class
 *
 * f(x|x_0, \gamma) = \gamma / [\pi (\gamma^2 + (x - x_0)^2)]
 *
 * X ~ Cauchy(x_0, \gamma)
 */
class RANDLIBSHARED_EXPORT CauchyRand : public ContinuousRand
{
    double x0, gamma;
    double gammaInv; /// 1 / gamma

public:
    CauchyRand(double location = 0, double scale = 1);
    std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return x0; }
    inline double getScale() const { return gamma; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double scale);
    static double standardVariate();

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;
};

#endif // CAUCHYRAND_H

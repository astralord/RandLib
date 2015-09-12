#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include "ContinuousRand.h"
#include "UniformRand.h"

/**
 * @brief The CauchyRand class
 */
class RANDLIBSHARED_EXPORT CauchyRand : public ContinuousRand
{
    double x0, gamma;
    double gammaInv; /// 1 / gamma

public:
    CauchyRand(double location = 0, double scale = 1);
    virtual std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return x0; }
    inline double getScale() const { return gamma; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double location, double scale);
    static double standardVariate();

    double E() const override;
    double Var() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double Entropy() const;
};

#endif // CAUCHYRAND_H

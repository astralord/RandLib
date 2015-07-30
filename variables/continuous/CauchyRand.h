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

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return x0; }
    inline double getScale() const { return gamma; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double variate() override;

    double E() const override { return NAN; }
    double Var() const override { return INFINITY; }

    inline double Entropy() const { return std::log(4 * gamma * M_PI); }
};

#endif // CAUCHYRAND_H

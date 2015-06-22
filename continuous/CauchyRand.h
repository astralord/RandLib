#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include <RandomVariable.h>
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT CauchyRand : public ContinuousRand
{
    double x0, gamma;
    double gammaInv; /// 1 / gamma
    UniformRand U;

public:
    CauchyRand(double location, double scale);

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return x0; }
    inline double getScale() const { return gamma; }

    virtual double pdf (double x) const override;
    virtual double cdf(double x) const override;
    virtual double value() override;

    inline double M() const override { return NAN; }
    inline double Var() const override { return INFINITY; }

    inline double Entropy() const { return std::log(4 * gamma * M_PI); }
};

#endif // CAUCHYRAND_H

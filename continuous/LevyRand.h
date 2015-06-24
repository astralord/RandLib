#ifndef LEVYRAND_H
#define LEVYRAND_H

#include <RandomVariable.h>
#include "NormalRand.h"

/**
 * @brief The LevyRand class
 * X ~ Levy(mu, c)
 * X ~ Stable(0.5, 1, mu, c)
 */
class RANDLIBSHARED_EXPORT LevyRand : public ContinuousRand
{
    double mu, c_2; /// c_2 = 0.5 * c
    double sqrtc_2pi; /// sqrtc_2pi = sqrt(0.5 * c / pi)
    NormalRand X;

public:
    LevyRand(double location = 0, double scale = 1);

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return c_2 + c_2; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double M() const override { return INFINITY; }
    double Var() const override { return INFINITY; }
};

#endif // LEVYRAND_H

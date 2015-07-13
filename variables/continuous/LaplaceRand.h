#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include "ContinuousRand.h"
#include "ExponentialRand.h"

/**
 * @brief The LaplaceRand class
 */
class RANDLIBSHARED_EXPORT LaplaceRand : public ContinuousRand
{
    double mu, b;
    double bInv; /// 1 / b
    ExponentialRand X;

public:
    LaplaceRand(double location = 0, double scale = 1);

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return b; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double E() const override { return mu; }
    double Var() const override { return 2 * b * b; }
};

#endif // LAPLACERAND_H

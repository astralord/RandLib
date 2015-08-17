#ifndef MAXWELLBOLTZMANNRAND_H
#define MAXWELLBOLTZMANNRAND_H

#include "ContinuousRand.h"
#include "ChiSquaredRand.h"

/**
 * @brief The MaxwellBoltzmannRand class
 */
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public ContinuousRand
{
    double a;

    ChiSquaredRand C;

public:
    MaxwellBoltzmannRand(double scale);
    virtual void setName() override;

    void setScale(double scale);
    inline double getScale() const { return a; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double variate() const override;

    double E() const override { return 2 * M_1_SQRTPI * M_SQRT2 * a; }
    double Var() const override { return a * a * (3 - 8.0 * M_1_PI); }
};

#endif // MAXWELLBOLTZMANNRAND_H

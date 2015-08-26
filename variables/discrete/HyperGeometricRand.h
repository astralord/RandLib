#ifndef HYPERGEOMETRICRAND_H
#define HYPERGEOMETRICRAND_H

#include "DiscreteRand.h"

class RANDLIBSHARED_EXPORT HyperGeometricRand : public DiscreteRand<int>
{
    unsigned N, K, n;
public:
    HyperGeometricRand();
    virtual void setName() override;

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return 1; }
    double Var() const override { return 1; }
};

#endif // HYPERGEOMETRICRAND_H

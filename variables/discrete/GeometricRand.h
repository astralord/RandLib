#ifndef GEOMETRICRAND_H
#define GEOMETRICRAND_H

#include "DiscreteRand.h"
#include "../continuous/ExponentialRand.h"

/**
 * @brief The GeometricRand class
 */
class RANDLIBSHARED_EXPORT GeometricRand : public DiscreteRand<int>
{
    double p;
    ExponentialRand Exp;

public:
    GeometricRand(double probability);
    inline void setProbability(double probability);
    inline double getProbability() { return p; }

    virtual double P(int k) const override;
    virtual double F(double x) const override;
    virtual int variate() override;

    double E() const override { return 1.0 / p; }
    double Var() const override { return (1 - p) / (p * p); }

    // TODO: add median and entropy
    static double constexpr Mode() { return 1; }
    inline double Skewness() { return (2 - p) / std::sqrt(1 - p); }
    inline double ExcessiveKurtosis() { return 1.0 / (p * (1 - p)) - 6; }
};

#endif // GEOMETRICRAND_H

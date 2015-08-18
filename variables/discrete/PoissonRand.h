#ifndef POISSONRAND_H
#define POISSONRAND_H

#include "DiscreteRand.h"
#include "../continuous/UniformRand.h"

/**
 * @brief The PoissonRand class
 */
class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteRand<int>
{
    double l;
    double exp_l; /// exp(-l)
public:
    PoissonRand(double rate);
    virtual void setName() override;

    void setRate(double rate);
    inline double getRate() const { return l; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return l; }
    double Var() const override { return l; }

    inline double Mode() const { return std::floor(l); }
    inline double Skewness() const { return 1.0 / std::sqrt(l); }
    inline double ExcessiveKurtosis() const { return 1.0 / l; }
};

#endif // POISSONRAND_H

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
    UniformRand U;
public:
    PoissonRand(double rate);

    void setRate(double rate);
    inline double getRate() const { return l; }

    virtual double P(int k) const override;
    virtual double F(double x) const override;
    virtual int value() override;

    double M() const override { return l; }
    double Var() const override { return l; }

    inline double Mode() const { return qFloor(l); } /// and also qCeil(l) - 1
    inline double Skewness() const { return 1.0 / std::sqrt(l); }
    inline double ExcessiveKurtosis() const { return 1.0 / l; }
};

#endif // POISSONRAND_H

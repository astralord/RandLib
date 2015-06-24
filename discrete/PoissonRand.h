#ifndef POISSONRAND_H
#define POISSONRAND_H

#include <RandomVariable.h>
#include <continuous/UniformRand.h>

class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteIntRand
{
    double l;
    double exp_l; // exp{-l}
    UniformRand U;
public:
    PoissonRand(double rate);

    void setRate(double rate);
     double getRate() const { return l; }

    virtual double P(int k) const override;
    virtual double F(double x) const override;
    virtual double value() override;

     double M() const override { return l; }
     double Var() const override { return l; }

     double Mode() const { return qFloor(l); } /// and also qCeil(l) - 1
     double Skewness() const { return 1.0 / std::sqrt(l); }
     double ExcessiveKurtosis() const { return 1.0 / l; }
};

#endif // POISSONRAND_H

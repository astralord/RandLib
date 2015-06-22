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
    inline double getRate() const { return l; }

    virtual double P(int k) const override;
    virtual double cdf(double x) const override;
    virtual double value() override;

    inline double M() const override { return l; }
    inline double Var() const override { return l; }

};

#endif // POISSONRAND_H

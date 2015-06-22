#ifndef POISSONRAND_H
#define POISSONRAND_H

#include <RandomVariable.h>
#include <continuous/UniformRand.h>

class RANDLIBSHARED_EXPORT PoissonRand : public DiscreteRand<int>
{
    double l;
    double exp_l; // exp{-l}
    UniformRand U;
public:
    PoissonRand(double rate);

    void setRate(double rate);
    double getRate() { return l; }

    virtual double P(int k);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return l; }
    inline double Var() { return l; }

};

#endif // POISSONRAND_H

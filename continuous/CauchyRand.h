#ifndef CAUCHYRAND_H
#define CAUCHYRAND_H

#include <RandomVariable.h>
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT CauchyRand : public ContinuousRandom
{
    double x0, gamma;
    double gammaInv; /// 1 / gamma
    UniformRand U;

public:
    CauchyRand(double location, double scale);

    void setLocation(double location);
    void setScale(double scale);
    double getLocation() { return x0; }
    double getScale() { return gamma; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();
    inline double M() { return NAN; }
    inline double Var() { return INFINITY; }

    inline double Entropy() { return qLn(4 * gamma * M_PI); }
};

#endif // CAUCHYRAND_H

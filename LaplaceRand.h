#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include <RandomVariable.h>
#include <ExponentialRand.h>

class RANDLIBSHARED_EXPORT LaplaceRand : public RandomVariable
{
    double mu, b;
    double bInv; /// 1 / b
    ExponentialRand X;

public:
    LaplaceRand(double location, double scale);

    void setLocation(double location);
    void setScale(double scale);
    double getLocation() { return mu; }
    double getScale() { return b; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return mu; }
    inline double Var() { return 2 * b * b; }
};

#endif // LAPLACERAND_H

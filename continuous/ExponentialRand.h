#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include <RandomVariable.h>
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT ExponentialRand : public ContinuousRandom
{
    double l, beta;

    //TODO: find a way to initialize them only once!!!
    /// Tables for ziggurat
    unsigned long ke[256];
    float we[256], fe[256];
    UniformRand U;
    double ziggurat();

public:
    ExponentialRand(double rate);

    void setRate(double rate);
    double getRate() { return l; }
    double getScale() { return beta; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();
    inline double M() { return beta; }
    inline double Var() { return beta * beta; }

    inline double Skewness() { return 2.0; }
    inline double ExcessiveKurtosis() { return 6.0; }

    inline double Entropy() { return 1 - qLn(l); }
};

#endif // EXPONENTIALRAND_H

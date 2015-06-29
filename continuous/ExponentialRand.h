#ifndef EXPONENTIALRAND_H
#define EXPONENTIALRAND_H

#include "ContinuousRand.h"
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT ExponentialRand : public ContinuousRand
{
    double l, beta;

    //TODO: find a way to initialize them only once!!!
    /// Tables for ziggurat
    unsigned long ke[256];
    float we[256], fe[256];
    UniformRand U;
    double ziggurat();

public:
    ExponentialRand(double rate = 1);

    void setRate(double rate);
    inline double getRate() const { return l; }
    inline double getScale() const { return beta; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double M() const override { return beta; }
    double Var() const override { return beta * beta; }

    inline double Skewness() const { return 2.0; }
    inline double ExcessiveKurtosis() const { return 6.0; }

    inline double Entropy() const { return 1 - std::log(l); }
};

#endif // EXPONENTIALRAND_H

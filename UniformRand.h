#ifndef UNIFORMRAND_H
#define UNIFORMRAND_H

#include <RandomVariable.h>

class RANDLIBSHARED_EXPORT UniformRand : public RandomVariable
{
    double a, b;
    double c; /// c = 1 / (b - a)
    void swapBoundaries();

public:
    UniformRand(double minValue, double maxValue);
    void setBoundaries(double minValue, double maxValue);

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return .5 * (b + a); }
    inline double Var() { return (b - a) * (b - a) / 12; }

    inline double Median() { return .5 * (b + a); }
    inline double Skewness() { return 0; }
    inline double ExcessKurtosis() { return -1.2; }

    inline double Entropy() { return (b == a) ? -INFINITY : qLn(b - a); }

};

#endif // UNIFORMRAND_H

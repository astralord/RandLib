#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include <RandomVariable.h>
#include <ChiSquaredRand.h>
#include <NormalRand.h>

class RANDLIBSHARED_EXPORT StudentTRand : public RandomVariable
{
    int v;
    NormalRand X;
    ChiSquaredRand Y;
public:
    StudentTRand(int degree);

    void setDegree(int degree);
    int getDegree() { return v; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return (v > 1) ? 0 : NAN; }
    inline double Var() {
        if (v > 2)
            return (double)v / (v - 2);
        return (v > 1) ? INFINITY : NAN;
    }
};

#endif // STUDENTTRAND_H

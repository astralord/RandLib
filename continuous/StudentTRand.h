#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include <RandomVariable.h>
#include "ChiSquaredRand.h"
#include "NormalRand.h"

class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousRand
{
    int v;
    NormalRand X;
    ChiSquaredRand Y;
    double pdfCoef;
public:
    StudentTRand(int degree);

    void setDegree(double degree);
    inline double getDegree() const { return v; }

    virtual double pdf (double x) const override;
    virtual double cdf(double x) const override;
    virtual double value() override;

    inline double M() const override { return (v > 1) ? 0 : NAN; }
    inline double Var() const override {
        if (v > 2)
            return (double)v / (v - 2);
        return (v > 1) ? INFINITY : NAN;
    }
};

#endif // STUDENTTRAND_H

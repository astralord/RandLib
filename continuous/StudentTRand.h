#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include <RandomVariable.h>
#include "ChiSquaredRand.h"
#include "NormalRand.h"

//TODO: degree should be int or double??
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousRand
{
    double v;
    NormalRand X;
    ChiSquaredRand Y;
    double pdfCoef;
public:
    StudentTRand(int degree);

    void setDegree(double degree);
    inline double getDegree() const { return v; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double M() const override { return (v > 1) ? 0 : NAN; }
    double Var() const override {
        if (v > 2)
            return v / (v - 2);
        return (v > 1) ? INFINITY : NAN;
    }
};

#endif // STUDENTTRAND_H

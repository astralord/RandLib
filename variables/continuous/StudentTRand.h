#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include "ContinuousRand.h"
#include "ChiSquaredRand.h"
#include "NormalRand.h"

/**
 * @brief The StudentTRand class
 */
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousRand
{
    int v;
    ChiSquaredRand Y;
    double pdfCoef;
public:
    StudentTRand(int degree);

    void setDegree(int degree);
    inline int getDegree() const { return v; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double variate() override;

    double E() const override { return (v > 1) ? 0 : NAN; }
    double Var() const override {
        if (v > 2)
            return v / (v - 2);
        return (v > 1) ? INFINITY : NAN;
    }
};

#endif // STUDENTTRAND_H

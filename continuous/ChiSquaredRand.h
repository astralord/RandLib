#ifndef CHISQUAREDRAND_H
#define CHISQUAREDRAND_H

#include <RandomVariable.h>
#include "NormalRand.h"

class RANDLIBSHARED_EXPORT ChiSquaredRand : public ContinuousRand
{
    int k;
    double pdfCoef, cdfCoef;
    NormalRand X;

public:
    ChiSquaredRand(int k);

    void setDegree(int degree);
    inline size_t getDegree() const { return k; }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    inline double M() const override { return k; }
    inline double Var() const override { return k + k; }

    inline double Mode() const { return std::max(k - 2, 0); }
    inline double Skewness() const { return std::sqrt(8.0 / k); }
    inline double ExcessKurtosis() const { return 12.0 / k; }
};

#endif // CHISQUAREDRAND_H

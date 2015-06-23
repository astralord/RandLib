#ifndef LOGNORMALRAND_H
#define LOGNORMALRAND_H

#include <RandomVariable.h>
#include "NormalRand.h"

class RANDLIBSHARED_EXPORT LogNormalRand : public ContinuousRand
{
    double expMu, expSigmaSq;
    NormalRand X;

public:
    LogNormalRand(double location, double scale);
    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return X.M(); }
    inline double getScale() const { return X.getSigma(); }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    inline double M() const override { return expMu * std::sqrt(expSigmaSq); }
    inline double Var() const override { return (expSigmaSq - 1) * expMu * expMu * expSigmaSq; }

    inline double Median() const { return expMu; }
    inline double Mode() const{ return expMu / expSigmaSq; }
    inline double Skewness() const { return (expSigmaSq + 2) * std::sqrt(expSigmaSq - 1); }
    inline double ExcessKurtosis() const;
};

#endif // LOGNORMALRAND_H

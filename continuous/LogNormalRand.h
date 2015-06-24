#ifndef LOGNORMALRAND_H
#define LOGNORMALRAND_H

#include <RandomVariable.h>
#include "NormalRand.h"

class RANDLIBSHARED_EXPORT LogNormalRand : public ContinuousRand
{
    double expMu, expVar;
    NormalRand X;

public:
    LogNormalRand(double location = 0, double scale = 1);
    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return X.M(); }
    inline double getScale() const { return X.getSigma(); }

    virtual double f(double x) const override;
    virtual double F(double x) const override;
    virtual double value() override;

    double M() const override { return expMu * std::sqrt(expVar); }
    double Var() const override { return (expVar - 1) * expMu * expMu * expVar; }

    inline double Median() const { return expMu; }
    inline double Mode() const{ return expMu / expVar; }
    inline double Skewness() const { return (expVar + 2) * std::sqrt(expVar - 1); }
    inline double ExcessKurtosis() const;

    bool fitToData(const QVector<double> &sample);
};

#endif // LOGNORMALRAND_H

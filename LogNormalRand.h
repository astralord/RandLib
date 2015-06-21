#ifndef LOGNORMALRAND_H
#define LOGNORMALRAND_H

#include <RandomVariable.h>
#include <NormalRand.h>

class RANDLIBSHARED_EXPORT LogNormalRand : public RandomVariable
{
    double expMu, expSigmaSq;
    NormalRand X;

public:
    LogNormalRand(double location, double scale);
    void setLocation(double location);
    void setScale(double scale);
    double getLocation() { return X.M(); }
    double getScale() { return X.getSigma(); }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return expMu * std::sqrt(expSigmaSq); }
    inline double Var() { return (expSigmaSq - 1) * expMu * expMu * expSigmaSq; }

    inline double Median() { return expMu; }
    inline double Mode() { return expMu / expSigmaSq; }
    inline double Skewness() { return (expSigmaSq + 2) * std::sqrt(expSigmaSq - 1); }
    inline double ExcessKurtosis();
};

#endif // LOGNORMALRAND_H

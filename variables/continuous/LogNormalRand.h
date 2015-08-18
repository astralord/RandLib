#ifndef LOGNORMALRAND_H
#define LOGNORMALRAND_H

#include "ContinuousRand.h"
#include "NormalRand.h"

/**
 * @brief The LogNormalRand class
 */
class RANDLIBSHARED_EXPORT LogNormalRand : public ContinuousRand
{
    double expMu, expVar;
    NormalRand X;

public:
    LogNormalRand(double location = 0, double scale = 1);
    virtual void setName() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return X.E(); }
    inline double getScale() const { return X.getSigma(); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return expMu * std::sqrt(expVar); }
    double Var() const override { return (expVar - 1) * expMu * expMu * expVar; }

    inline double Median() const { return expMu; }
    inline double Mode() const{ return expMu / expVar; }
    inline double Skewness() const { return (expVar + 2) * std::sqrt(expVar - 1); }
    inline double ExcessKurtosis() const;

    bool fitToData(const QVector<double> &sample);
};

#endif // LOGNORMALRAND_H

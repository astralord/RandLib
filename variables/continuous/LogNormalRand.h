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
    virtual std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return X.E(); }
    inline double getScale() const { return X.getSigma(); }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override;
    double Var() const override;

    double Quantile(double p);

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool fitToData(const QVector<double> &sample);
};

#endif // LOGNORMALRAND_H

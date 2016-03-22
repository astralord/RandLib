#ifndef LOGISTICRAND_H
#define LOGISTICRAND_H

#include "ContinuousDistribution.h"
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT LogisticRand : public ContinuousDistribution
{
    double mu, s;

public:
    LogisticRand(double location = 0, double scale = 1);
    std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return s; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const { return mu; }
    double Variance() const;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    /// Method of moments
    bool fitLocation_MM(const QVector<double> &sample);
    bool fitScale_MM(const QVector<double> &sample);
    bool fit_MM(const QVector<double> &sample);
};

#endif // LOGISTICRAND_H

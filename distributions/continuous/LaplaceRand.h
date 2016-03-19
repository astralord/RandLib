#ifndef LAPLACERAND_H
#define LAPLACERAND_H

#include "ExponentialRand.h"

/**
 * @brief The LaplaceRand class
 */
class RANDLIBSHARED_EXPORT LaplaceRand : public ContinuousDistribution
{
    double mu, b;
    double bInv; /// 1 / b

public:
    LaplaceRand(double location = 0, double scale = 1);
    std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return b; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double location, double scale);

    double Mean() const override;
    double Variance() const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    std::complex<double> CF(double t) const override;

    /// Maximum likelihood estimation
    bool fitLocation_MLE(const QVector<double> &sample);
    bool fitScale_MLE(const QVector<double> &sample);
    bool fit_MLE(const QVector<double> &sample);
    
    /// Method of moments
    bool fitLocation_MM(const QVector<double> &sample);
    bool fitScale_MM(const QVector<double> &sample);
    bool fit_MM(const QVector<double> &sample);
};

#endif // LAPLACERAND_H

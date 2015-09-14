#ifndef RAYLEIGHRAND_H
#define RAYLEIGHRAND_H

#include "ContinuousRand.h"
#include "ExponentialRand.h"

/**
 * @brief The RayleighRand class
 */
class RANDLIBSHARED_EXPORT RayleighRand : public ContinuousRand
{
    double sigma;
    double sigmaSqInv; // 1.0 / sigma^2

public:
    explicit RayleighRand(double scale = 1);
    virtual std::string name() override;

    void setScale(double scale);
    inline double getScale() const { return sigma; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool fitToData(const QVector<double> &sample);
};
#endif // RAYLEIGHRAND_H

#ifndef LOGCAUCHYRAND_H
#define LOGCAUCHYRAND_H

#include "ContinuousDistribution.h"
#include "CauchyRand.h"

/**
 * @brief The LogCauchyRand class
 */
class RANDLIBSHARED_EXPORT LogCauchyRand : public ContinuousDistribution
{
    double expMu;
    CauchyRand X;

public:
    LogCauchyRand(double location = 0, double scale = 1);
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return SEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return X.getLocation(); }
    inline double getScale() const { return X.getScale(); }

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
};

#endif // LOGCAUCHYRAND_H

#ifndef LOGCAUCHYRAND_H
#define LOGCAUCHYRAND_H

#include "ContinuousDistribution.h"
#include "CauchyRand.h"

/**
 * @brief The LogCauchyRand class
 */
class RANDLIBSHARED_EXPORT LogCauchyRand : public ContinuousDistribution
{
    CauchyRand X;

public:
    LogCauchyRand(double location = 0, double scale = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return X.GetLocation(); }
    inline double GetScale() const { return X.GetScale(); }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

public:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};

#endif // LOGCAUCHYRAND_H

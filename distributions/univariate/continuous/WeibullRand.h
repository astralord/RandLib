#ifndef WEIBULLRAND_H
#define WEIBULLRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The WeibullRand class
 */
class RANDLIBSHARED_EXPORT WeibullRand : public ContinuousDistribution
{
    double lambda, k;
    double kInv;

public:
    WeibullRand(double scale = 1, double shape = 1);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void setParameters(double scale, double shape);
    inline double getScale() const { return lambda; }
    inline double getShape() const { return k; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double QuantileImpl(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;
};

#endif // WEIBULLRAND_H

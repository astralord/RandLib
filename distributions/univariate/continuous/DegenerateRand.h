#ifndef DEGENERATERAND_H
#define DEGENERATERAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The DegenerateRand class
 * Degenerate distribution
 * X ~ \delta(a)
 *
 * f(x|a) = \delta(a)
 */
class RANDLIBSHARED_EXPORT DegenerateRand : public ContinuousDistribution
{
    double a;

public:
    explicit DegenerateRand(double value = 0);
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return a; }
    double MaxValue() const override { return a; }

    void setValue(double value);
    inline double getValue() const { return a; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;
};

#endif // DEGENERATERAND_H

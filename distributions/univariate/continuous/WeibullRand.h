#ifndef WEIBULLRAND_H
#define WEIBULLRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The WeibullRand class
 * Weibull distribution
 */
class RANDLIBSHARED_EXPORT WeibullRand : public ContinuousDistribution
{
    double lambda, k;
    double kInv;

public:
    WeibullRand(double scale = 1, double shape = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double scale, double shape);
    inline double GetScale() const { return lambda; }
    inline double GetShape() const { return k; }

    double f(double x) const override;
    double logf(double x) const override;
    double F(double x) const override;
    double S(double x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    double Entropy() const;
};

#endif // WEIBULLRAND_H

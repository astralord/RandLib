#ifndef SECHRAND_H
#define SECHRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The SechRand class
 */
class RANDLIBSHARED_EXPORT SechRand : public ContinuousDistribution
{
public:
    SechRand();
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return INFINITE_T; }
    double MinValue() const override { return INFINITY; }
    double MaxValue() const override { return INFINITY; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double Kurtosis() const override;

    double Entropy() const;
};

#endif // SECHRAND_H

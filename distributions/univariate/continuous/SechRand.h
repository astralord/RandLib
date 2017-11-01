#ifndef SECHRAND_H
#define SECHRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The SechRand class <BR>
 * Hyperbolic secant distribution
 *
 * Notation: X ~ Sech
 */
class RANDLIBSHARED_EXPORT SechRand : public ContinuousDistribution
{
public:
    SechRand();

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const;
    double quantileImpl1m(double p) const;

    std::complex<double> CFImpl(double t) const override;
public:
    double Entropy() const;
};

#endif // SECHRAND_H

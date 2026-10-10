#ifndef SECHRAND_H
#define SECHRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The SechRand class <BR>
 * Hyperbolic secant distribution
 *
 * Notation: X ~ Sech
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT SechRand : public ContinuousDistribution<RealType>
{
public:
    SechRand();

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    RealType MinValue() const override { return -INFINITY; }
    RealType MaxValue() const override { return INFINITY; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    RealType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;
public:
    long double Entropy() const;
};

#endif // SECHRAND_H

#ifndef KOLMOGOROVSMIRNOVRAND_H
#define KOLMOGOROVSMIRNOVRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The KolmogorovSmirnovRand class <BR>
 * Kolmogorov-Smirnov distribution
 *
 * Notation: X ~ KS
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT KolmogorovSmirnovRand : public ContinuousDistribution<RealType>
{
public:
    KolmogorovSmirnovRand();
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

private:
    static double L(RealType x);
    static double K(RealType x);
public:
    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    double logF(const RealType & x) const;
    double logS(const RealType & x) const;

private:
    RealType truncatedGammaVariate() const;
    RealType variateForTheLeftMostInterval() const;
    RealType variateForTheRightMostInterval() const;
public:
    RealType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Mode() const override;
    RealType Median() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
};

#endif // KOLMOGOROVSMIRNOVRAND_H

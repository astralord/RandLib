#ifndef KOLMOGOROVSMIRNOVRAND_H
#define KOLMOGOROVSMIRNOVRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The KolmogorovSmirnovRand class <BR>
 * Kolmogorov-Smirnov distribution
 *
 * Notation: X ~ KS
 */
class RANDLIBSHARED_EXPORT KolmogorovSmirnovRand : public ContinuousDistribution
{
public:
    KolmogorovSmirnovRand();
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

private:
    static double L(double x);
    static double K(double x);
public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double logF(const double & x) const;
    double logS(const double & x) const;

private:
    double truncatedGammaVariate() const;
    double variateForTheLeftMostInterval() const;
    double variateForTheRightMostInterval() const;
public:
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Median() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};

#endif // KOLMOGOROVSMIRNOVRAND_H

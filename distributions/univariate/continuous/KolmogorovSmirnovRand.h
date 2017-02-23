#ifndef KOLMOGOROVSMIRNOVRAND_H
#define KOLMOGOROVSMIRNOVRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The KolmogorovSmirnovRand class
 */
class RANDLIBSHARED_EXPORT KolmogorovSmirnovRand : public ContinuousDistribution
{
public:
    KolmogorovSmirnovRand();
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

private:
    static double PDF(double x);
    static double L(double x);
    static double K(double x);
    static double CDF(double x);
    static double CDFCompl(double x);
public:
    double f(double x) const override;
    double logf(double x) const override;
    double F(double x) const override;
    double S(double x) const override;

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

public:
    static double Quantile(double p);
    static double Quantile1m(double p);
private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
};

#endif // KOLMOGOROVSMIRNOVRAND_H

#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include "GammaRand.h"

/**
 * @brief The StudentTRand class
 * Student's t-distribution
 *
 * Notation: X ~ t(nu)
 */
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousDistribution
{
    double nu;
    double mu, sigma;
    GammaRand Y;
    double pdfCoef;
    double nup1Half; /// 0.5 * (nu + 1)
    double betaInv; /// 1.0 / beta(0.5 * nu, 0.5)

public:
    explicit StudentTRand(double degree, double location = 0.0, double scale = 1.0);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetDegree(double degree);
    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetDegree() const { return nu; }
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return sigma; }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

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
};

#endif // STUDENTTRAND_H

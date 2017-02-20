#ifndef INVERSEGAUSSIANRAND_H
#define INVERSEGAUSSIANRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The InverseGaussianRand class
 * Inverse Gaussian (Wald) distribution
 *
 * Notation: X ~ IG(μ, λ)
 */
class RANDLIBSHARED_EXPORT InverseGaussianRand : public ContinuousDistribution
{
    double mu, lambda;

    double pdfCoef; /// (λ/(2π))^(1/2)
    double cdfCoef; /// exp(2λ/μ)
public:
    InverseGaussianRand(double mean = 1, double shape = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double mean, double shape);
    inline double GetMean() const { return mu; }
    inline double GetShape() const { return lambda; }

    double f(double x) const override;
    double F(double x) const override;
    double S(double x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};


#endif // INVERSEGAUSSIANRAND_H

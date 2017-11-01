#ifndef INVERSEGAUSSIANRAND_H
#define INVERSEGAUSSIANRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The InverseGaussianRand class <BR>
 * Inverse Gaussian (Wald) distribution
 *
 * Notation: X ~ IG(μ, λ)
 */
class RANDLIBSHARED_EXPORT InverseGaussianRand : public ContinuousDistribution
{
    double mu = 1; ///< mean μ
    double lambda = 1; ///< shape λ
    double pdfCoef = M_SQRT1_2 / M_SQRTPI; ///< (λ/(2π))^(1/2)
    double cdfCoef = M_E * M_E; ///< exp(2λ/μ)
public:
    InverseGaussianRand(double mean = 1, double shape = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double mean, double shape);
    inline double GetMean() const { return mu; }
    inline double GetShape() const { return lambda; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
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

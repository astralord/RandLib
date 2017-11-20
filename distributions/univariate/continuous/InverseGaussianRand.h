#ifndef INVERSEGAUSSIANRAND_H
#define INVERSEGAUSSIANRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The InverseGaussianRand class <BR>
 * Inverse Gaussian (Wald) distribution
 *
 * Notation: X ~ IG(μ, λ)
 */
template < typename RealType = long double >
class RANDLIBSHARED_EXPORT InverseGaussianRand : public ContinuousDistribution<RealType>
{
    double mu = 1; ///< mean μ
    double lambda = 1; ///< shape λ
    double pdfCoef = M_SQRT1_2 / M_SQRTPI; ///< (λ/(2π))^(1/2)
    double cdfCoef = M_E * M_E; ///< exp(2λ/μ)
public:
    InverseGaussianRand(double mean = 1, double shape = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double mean, double shape);
    inline double GetMean() const { return mu; }
    inline double GetShape() const { return lambda; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};


#endif // INVERSEGAUSSIANRAND_H

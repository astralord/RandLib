#ifndef STUDENTTRAND_H
#define STUDENTTRAND_H

#include "NakagamiRand.h"

/**
 * @brief The StudentTRand class <BR>
 * Student's t-distribution
 *
 * Notation: X ~ t(ν, μ, σ)
 * If X ~ t(1, μ, σ), then X ~ Cauchy(μ, σ)
 * X -> Normal(μ, σ) for t -> ∞
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT StudentTRand : public ContinuousDistribution<RealType>
{
    double nu = 1; ///< degree ν
    double mu = 0; ///< location μ
    double sigma = 1; ///< scale σ
    double logSigma = 0; ///< log(σ)
    NakagamiRand<RealType> Y{};
    double pdfCoef = -M_LNPI; ///< coefficient for faster pdf calculation
    double nup1Half = 1; ///< 0.5 * (ν + 1)
    double logBetaFun = M_LNPI; ///< log(B(0.5 * ν, 0.5))

public:
    explicit StudentTRand(double degree = 1.0, double location = 0.0, double scale = 1.0);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    RealType MinValue() const override { return -INFINITY; }
    RealType MaxValue() const override { return INFINITY; }

    void SetDegree(double degree);
    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetDegree() const { return nu; }
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return sigma; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

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
};

#endif // STUDENTTRAND_H

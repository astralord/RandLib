#ifndef EXPONENTIALLYMODIFIEDGAUSSIANRAND_H
#define EXPONENTIALLYMODIFIEDGAUSSIANRAND_H

#include "ContinuousDistribution.h"
#include "NormalRand.h"
#include "ExponentialRand.h"

/**
 * @brief The ExponentiallyModifiedGaussianRand class <BR>
 * Exponentially-modified Gaussian distribution
 *
 * Notation: X ~ EMG(μ, σ, β)
 *
 * Related distributions: <BR>
 * X = Y + Z, where Y ~ Normal(μ, σ) and Z ~ Exp(β)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT ExponentiallyModifiedGaussianRand : public ContinuousDistribution<RealType>
{
    NormalRand<RealType> X{};
    ExponentialRand<RealType> Y{};

    double a = 1; ///< μ + βσ^2
    double b = M_SQRT1_2; ///< 1 / (√2 * σ)
    double c = 0.5; ///< μ + βσ^2 / 2
    double v = 1; /// βσ

public:
    explicit ExponentiallyModifiedGaussianRand(double location = 0, double variance = 1, double rate = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    RealType MinValue() const override { return -INFINITY; }
    RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double location, double variance, double rate);
    inline double GetLocation() const { return X.GetLocation(); }
    inline double GetScale() const { return X.GetScale(); }
    inline double GetRate() const { return Y.GetRate(); }

private:
    DoublePair faux(const RealType &x) const;

public:
    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;
    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // EXPONENTIALLYMODIFIEDGAUSSIANRAND_H

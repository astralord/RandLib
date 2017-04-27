#ifndef EXPONENTIALLYMODIFIEDGAUSSIANRAND_H
#define EXPONENTIALLYMODIFIEDGAUSSIANRAND_H

#include "ContinuousDistribution.h"
#include "NormalRand.h"
#include "ExponentialRand.h"

/**
 * @brief The ExponentiallyModifiedGaussianRand class
 * Exponentially-modified Gaussian distribution
 *
 * Notation: X ~ EMG(μ, σ, β)
 *
 * Related distributions:
 * X = Y + Z, where Y ~ Normal(μ, σ) and Z ~ Exp(β)
 */
class RANDLIBSHARED_EXPORT ExponentiallyModifiedGaussianRand : public ContinuousDistribution
{
    NormalRand X;
    ExponentialRand Y;

    /// auxiliary variables
    double a, b, c, v;

public:
    explicit ExponentiallyModifiedGaussianRand(double location = 0, double variance = 1, double rate = 1);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double location, double variance, double rate);
    inline double GetLocation() const { return X.GetLocation(); }
    inline double GetScale() const { return X.GetScale(); }
    inline double GetRate() const { return Y.GetRate(); }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;
    static double StandardVariate();

    double Mean() const override;
    double Variance() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // EXPONENTIALLYMODIFIEDGAUSSIANRAND_H

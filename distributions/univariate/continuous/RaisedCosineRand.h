#ifndef RAISEDCOSINERAND_H
#define RAISEDCOSINERAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The RaisedCosineDistribution class <BR>
 * Abstract class for Raised-cosine distribution
 *
 * Notation: X ~ Raised-cosine(μ, s)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT RaisedCosineDistribution : public ContinuousDistribution<RealType>
{
    double mu = 0; ///< location μ
    double s = M_PI; ///< scale
    double s_pi = 1; ///< s / π
    double log2S = M_LN2 + M_LNPI; ///< log(2s)

protected:
    RaisedCosineDistribution(double location, double scale);

public:
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    RealType MinValue() const override { return mu - s; }
    RealType MaxValue() const override { return mu + s; }

protected:
    void SetLocation(double location);
    void SetScale(double scale);

public:
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return s; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

    static RealType StandardVariate(RandGenerator &randGenerator = ProbabilityDistribution<RealType>::staticRandGenerator);
    RealType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

/**
 * @brief The RaisedCosineRand class <BR>
 * Raised-cosine distribution
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT RaisedCosineRand : public RaisedCosineDistribution<RealType>
{
public:
    RaisedCosineRand(double location = 0, double scale = M_PI) : RaisedCosineDistribution<RealType>(location, scale) {}
    String Name() const override;

    using RaisedCosineDistribution<RealType>::SetLocation;
    using RaisedCosineDistribution<RealType>::SetScale;
};

/**
 * @brief The RaabGreenRand class <BR>
 * Raab-Green distribution
 *
 * Notation: X ~ Raab-Green()
 *
 * Related distributions:
 * X ~ Raised-cosine(0.0, π)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT RaabGreenRand : public RaisedCosineDistribution<RealType>
{
public:
    RaabGreenRand() : RaisedCosineDistribution<RealType>(0.0, M_PI) {}
    String Name() const override;
};


#endif // RAISEDCOSINERAND_H

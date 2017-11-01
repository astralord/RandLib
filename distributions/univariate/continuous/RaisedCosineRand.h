#ifndef RAISEDCOSINERAND_H
#define RAISEDCOSINERAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The RaisedCosineDistribution class <BR>
 * Abstract class for Raised-cosine distribution
 *
 * Notation: X ~ Raised-cosine(μ, s)
 */
class RANDLIBSHARED_EXPORT RaisedCosineDistribution : public ContinuousDistribution
{
    double mu = 0; ///< location μ
    double s = M_PI; ///< scale
    double s_pi = 1; ///< s / π
    double log2S = M_LN2 + M_LNPI; ///< log(2s)

protected:
    RaisedCosineDistribution(double location, double scale);

public:
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return mu - s; }
    double MaxValue() const override { return mu + s; }

protected:
    void SetLocation(double location);
    void SetScale(double scale);

public:
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return s; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;

    static double StandardVariate();
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

/**
 * @brief The RaisedCosineRand class <BR>
 * Raised-cosine distribution
 */
class RANDLIBSHARED_EXPORT RaisedCosineRand : public RaisedCosineDistribution
{
public:
    RaisedCosineRand(double location, double scale) : RaisedCosineDistribution(location, scale) {}
    String Name() const override;

    using RaisedCosineDistribution::SetLocation;
    using RaisedCosineDistribution::SetScale;
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
class RANDLIBSHARED_EXPORT RaabGreenRand : public RaisedCosineDistribution
{
public:
    RaabGreenRand() : RaisedCosineDistribution(0.0, M_PI) {}
    String Name() const override;
};


#endif // RAISEDCOSINERAND_H

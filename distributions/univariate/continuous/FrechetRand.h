#ifndef FRECHETRAND_H
#define FRECHETRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The FrechetRand class <BR>
 * Frechet distribution
 *
 * Notation: X ~ Frechet(α, s, m)
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT FrechetRand : public ContinuousDistribution<RealType>
{
    double alpha = 1; ///< shape α
    double s = 1; ///< scale
    double m = 0; ///< location
    double alphaInv = 1; ///< 1/α
    double pdfCoef = 0; ///< log(α/s)

public:
    FrechetRand(double shape = 1, double scale = 1, double location = 0);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return m; }
    RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double scale, double location);
    inline double GetShape() const { return alpha; }
    inline double GetScale() const { return s; }
    inline double GetLocation() const { return m; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;
    RealType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

public:
    long double Entropy() const override;
};

#endif // FRECHETRAND_H

#ifndef FRECHETRAND_H
#define FRECHETRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The FrechetRand class <BR>
 * Frechet distribution
 *
 * Notation: X ~ Frechet(α, s, m)
 */
class RANDLIBSHARED_EXPORT FrechetRand : public ContinuousDistribution
{
    double alpha = 1; ///< shape α
    double s = 1; ///< scale
    double m = 0; ///< location
    double alphaInv = 1; ///< 1/α
    double pdfCoef = 0; ///< log(α/s)

public:
    FrechetRand(double shape, double scale, double location);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return m; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double scale, double location);
    inline double GetShape() const { return alpha; }
    inline double GetScale() const { return s; }
    inline double GetLocation() const { return m; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double S(const double & x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    double Entropy() const;
};

#endif // FRECHETRAND_H

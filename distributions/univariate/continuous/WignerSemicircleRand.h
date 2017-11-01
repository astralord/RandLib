#ifndef WIGNERSEMICIRCLERAND_H
#define WIGNERSEMICIRCLERAND_H

#include "BetaRand.h"

/**
 * @brief The WignerSemicircleRand class <BR>
 * Wigner-Semicircle distribution
 *
 * Notation: X ~ Wigner-Sc(R)
 *
 * Related distributions:
 * If Y ~ Beta(1.5, 1.5), then R * (2Y - 1) ~ Wigner-Sc(R)
 */
class RANDLIBSHARED_EXPORT WignerSemicircleRand : public ContinuousDistribution
{
    double R = 1; ///< radius
    double RSq = 1; ///< R^2
    double logRSq = 0; /// log(R^2)
    BetaRand X{1.5, 1.5};
    
public:
    explicit WignerSemicircleRand(double radius);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return -R; }
    double MaxValue() const override { return R; }

    void SetRadius(double radius);
    inline double GetRadius() const { return R; }
    
public:
    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    double Entropy() const;
};

#endif // WIGNERSEMICIRCLERAND_H

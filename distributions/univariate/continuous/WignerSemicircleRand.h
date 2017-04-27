#ifndef WIGNERSEMICIRCLERAND_H
#define WIGNERSEMICIRCLERAND_H

#include "BetaRand.h"

/**
 * @brief The WignerSemicircleRand class
 * Wigner-Semicircle distribution
 *
 * Notation: X ~ Wigner-Sc(R)
 *
 * Related distributions:
 * If Y ~ Beta(1.5, 1.5), then R * (2Y - 1) ~ Wigner-Sc(R)
 */
class RANDLIBSHARED_EXPORT WignerSemicircleRand : public ContinuousDistribution
{
    double R, RSq;
    double logRSq; /// log(R^2)
    BetaRand X; /// for generator
    
public:
    explicit WignerSemicircleRand(double radius);

    std::string Name() const override;
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

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
template < typename RealType = double >
class RANDLIBSHARED_EXPORT WignerSemicircleRand : public ContinuousDistribution<RealType>
{
    RealType R = 1; ///< radius
    double RSq = 1; ///< R^2
    double logRSq = 0; /// log(R^2)
    BetaRand<RealType> X{1.5, 1.5};
    
public:
    explicit WignerSemicircleRand(double radius);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    RealType MinValue() const override { return -R; }
    RealType MaxValue() const override { return R; }

    void SetRadius(double radius);
    inline double GetRadius() const { return R; }
    
public:
    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    RealType Variate() const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;
    
    long double Entropy() const override;
};

#endif // WIGNERSEMICIRCLERAND_H

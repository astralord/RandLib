#ifndef STABLECOUNTRAND_H
#define STABLECOUNTRAND_H

#include "StableRand.h"

/**
 * @brief The StableCountRand class <BR>
 * Stable count distribution
 *
 * Notation: X ~ S-c(α, θ, ν)
 *
 * Related distributions: <BR>
 * If X ~ S(α, 1, γ, 0), where then γ = \cos(\pi α / 2) ** (1 / α), <BR>
 * then ...
 * TODO: verify this relation
 * If X ~ S-c(1/2, 1, θ), then X ~ Γ(3/2, 1/(4θ))
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT StableCountRand : public ContinuousDistribution<RealType>
{
private:
    StableRand<RealType> X{};

    double alpha = 0.5; ///< characteristic exponent α
    double theta = 1; ///< scale θ
    double nu = 0; ///< location ν

    double logAlpha = -M_LN2; ///< log(α)
    double gammaAlphaInv = 1; ///< Γ(1/α)
    double lgammaAlphaInv = 0; ///< log(Γ(1/α))
    double lgamma2AlphaInv = M_LN3 + M_LN2; ///< log(Γ(2/α))
    double logTheta = 0; ///< log(θ)

public:
    StableCountRand(double exponent = 0.5, double scale = 1, double location = 0);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    RealType MinValue() const override { return nu; }
    RealType MaxValue() const override { return INFINITY; }

    void SetExponent(double exponent);
    void SetLocation(double location);
    void SetScale(double scale);

    /**
     * @fn GetExponent
     * @return characteristic exponent α
     */
    inline double GetExponent() const { return alpha; }
    /**
     * @fn GetScale
     * @return scale parameter θ
     */
    inline double GetScale() const { return theta; }
    /**
     * @fn GetLocation
     * @return location parameter ν
     */
    inline double GetLocation() const { return nu; }

private:
    double logf_adj(const RealType &x) const;
    double f_adj(const RealType &x) const;
    double F_adj(const RealType &x) const;
    double S_adj(const RealType &x) const;
    RealType logsumexp(std::vector<RealType> &sample) const;

public:
    double logf(const RealType & x) const override;
    double f(const RealType &x) const override;
    double F(const RealType & x) const override;
    double S(const RealType &x) const override;
    RealType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
};

#endif // STABLECOUNTRAND_H

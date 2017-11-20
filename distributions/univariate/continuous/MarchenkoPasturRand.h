#ifndef MARCHENKOPASTURRAND_H
#define MARCHENKOPASTURRAND_H

#include "ContinuousDistribution.h"
#include "BetaRand.h"
/**
 * @brief The MarchenkoPasturRand class <BR>
 * Marchenko-Pastur distribution
 *
 * Notation: X ~ Marchenko-Pastur(λ, σ)
 */
template < typename RealType = long double>
class RANDLIBSHARED_EXPORT MarchenkoPasturRand : public ContinuousDistribution<RealType>
{
    double lambda = 1; ///< ratio index λ
    double sigmaSq = 1; ///< scale parameter σ^2
    double a = 0; ///< minimal value (apart from 0 if λ > 1)
    double b = 4; ///< maximal value
    double logLambda = 0; /// < log(λ)

    BetaRand<RealType> BetaRV{0.5, 1.5, 0, 4}; ///< beta-distributed rv for generator
    double M = 1.0; ///< rejection constant

public:
    MarchenkoPasturRand(double ratio, double scale);
    String Name() const override;

    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    RealType MinValue() const override { return (lambda < 1) ? sigmaSq * a : 0; }
    RealType MaxValue() const override { return sigmaSq * b; }

    void SetParameters(double ratio, double scale);
    double GetRatio() const { return lambda; }
    double GetScale() const { return sigmaSq; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;

private:
    /**
     * @fn ccdfForLargeRatio
     * @param x
     * @return S(x) for λ > 1
     */
    double ccdfForLargeRatio(const RealType & x) const;
    /**
     * @fn cdfForSmallRatio
     * @param x
     * @return F(x) for 0 < λ <= 1
     */
    double cdfForSmallRatio(const RealType & x) const;

public:
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

private:
    enum GENERATOR_ID {
        TINY_RATIO,
        SMALL_RATIO,
        LARGE_RATIO,
        HUGE_RATIO
    };

    GENERATOR_ID getIdOfUsedGenerator() const
    {
        if (lambda < 0.3)
            return TINY_RATIO;
        if (lambda <= 1.0)
            return SMALL_RATIO;
        return (lambda > 3.3) ? HUGE_RATIO : LARGE_RATIO;
    }

    RealType variateForTinyRatio() const;
    RealType variateForSmallRatio() const;
    RealType variateForLargeRatio() const;
    RealType variateForHugeRatio() const;

public:
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

private:
    long double Moment(int n) const;

public:
    long double Mean() const override;
    long double Variance() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

#endif // MARCHENKOPASTURRAND_H

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
class RANDLIBSHARED_EXPORT MarchenkoPasturRand : public ContinuousDistribution
{
    double lambda = 1, sigmaSq = 1;
    double a = 0, b = 4;
    BetaRand BetaRV{0.5, 1.5, 0, 4};
    double M = 1.0; ///< rejection constant

public:
    MarchenkoPasturRand(double ratio, double scale);
    String Name() const override;
    void SetParameters(double ratio, double scale);
    double GetRatio() const { return lambda; }
    double GetScale() const { return sigmaSq; }

    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return (lambda < 1) ? sigmaSq * a : 0; }
    double MaxValue() const override { return sigmaSq * b; }

    double f(const double & x) const override;
    double logf(const double & x) const override;

private:
    /**
     * @fn ccdfForLargeRatio
     * @param x
     * @return S(x) for λ > 1
     */
    double ccdfForLargeRatio(const double & x) const;
    /**
     * @fn cdfForSmallRatio
     * @param x
     * @return F(x) for 0 < λ <= 1
     */
    double cdfForSmallRatio(const double & x) const;

public:
    double F(const double & x) const override;
    double S(const double & x) const override;

private:
    enum GENERATOR_ID {
        TINY_RATIO,
        SMALL_RATIO,
        LARGE_RATIO,
        HUGE_RATIO
    };

    GENERATOR_ID getIdOfUsedGenerator() const;

    double variateForTinyRatio() const;
    double variateForSmallRatio() const;
    double variateForLargeRatio() const;
    double variateForHugeRatio() const;

public:
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

private:
    double Moment(int n) const;

public:
    double Mean() const override;
    double Variance() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

#endif // MARCHENKOPASTURRAND_H

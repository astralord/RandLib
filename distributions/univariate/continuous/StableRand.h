#ifndef STABLERAND_H
#define STABLERAND_H

#include <functional>
#include "LimitingDistribution.h"

/**
 * @brief The StableRand class
 * Stable distribution
 *
 * X ~ Stable(α, β, σ, μ)
 *
 * Related distributions:
 * If X ~ Normal(μ, σ), then X ~ Stable(2, 0, σ, μ)
 * If X ~ Cauchy(μ, σ), then X ~ Stable(1, 0, σ / 2^(1/2), μ)
 * If +/-X ~ Levy(μ, σ), then X ~ Stable(0.5, +/-1, σ, μ)
 */
class RANDLIBSHARED_EXPORT StableRand : public LimitingDistribution
{
    double xi, S, zeta; /// coefficients for common α
    double alpham1Inv, alpha_alpham1; /// 1 / (α - 1) and α / (α - 1)

    static constexpr double BIG_NUMBER = 1e9; /// aka infinity for pdf and cdf calculations
    static constexpr double ALMOST_TWO = 1.99999;

    enum DISTRIBUTION_ID {
        NORMAL, /// α = 2
        LEVY, /// α = 0.5, |β| = 1
        CAUCHY, /// α = 1, β = 0
        UNITY_EXPONENT, /// α = 1, β != 0
        COMMON /// the rest
    };

    DISTRIBUTION_ID distributionId;

protected:
    double pdfCoef;
    double pdfXLimit; /// boundary k such that for |x| > k we can use tail approximation

public:
    StableRand(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~StableRand() {}

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override {
        if (alpha < 1) {
            if (beta == 1)
                return RIGHTSEMIFINITE_T;
            if (beta == -1)
                return LEFTSEMIFINITE_T;
        }
        return INFINITE_T;
    }
    double MinValue() const override { return (alpha < 1 && beta == 1) ? mu : -INFINITY; }
    double MaxValue() const override { return (alpha < 1 && beta == -1) ? mu : INFINITY; }

    void SetParameters(double exponent, double skewness);
    void SetScale(double scale);

    /// Probability distribution functions
protected:
    double pdfNormal(double x) const;
    double pdfCauchy(double x) const;
    double pdfLevy(double x) const;
private:
    static double fastpdfExponentiation(double u);

    /// functions for pdf calculation for α = 1
    double limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const;
    double integrandAuxForUnityExponent(double theta, double xAdj) const;
    double integrandForUnityExponent(double theta, double xAdj) const;
    double pdfForUnityExponent(double x) const;

    /// functions for pdf calculation for α != 1
    DoublePair seriesZeroParams;
    double pdfAtZero() const; /// f(0)
    double pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const;
    double pdfSeriesExpansionAtInf(double logX, double xiAdj, int k) const;
    double pdfTaylorExpansionTailNearCauchy(double x) const; /// ~ f(x, α) - f(x, 1)
    double limitCaseForIntegrandAuxForCommonExponent(double theta, double xiAdj) const;
    double integrandAuxForCommonExponent(double theta, double xAdj, double xiAdj) const;
    double integrandForCommonExponent(double theta, double xAdj, double xiAdj) const;
    double pdfForCommonExponent(double x) const;
public:    
    double f(double x) const override;

    /// Cumulative distribution functions
protected:
    double cdfNormal(double x) const;
    double cdfCauchy(double x) const;
    double cdfLevy(double x) const;
private:
    static double fastcdfExponentiation(double u);

    double cdfForUnityExponent(double x) const;
    double cdfIntegralRepresentation(double absXSt, double xiAdj) const;
    double cdfForCommonExponent(double x) const;
public:
    double F(double x) const override;

    /// Variates
private:
    double variateForUnityExponent() const;
    double variateForCommonExponent() const;
public:
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

public:
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;
};


/**
 * @brief The HoltsmarkRand class
 *
 * X ~ Holtsmark(σ, μ)
 *
 * Related distributions:
 * X ~ Stable(1.5, 0, σ, μ)
 */
class RANDLIBSHARED_EXPORT HoltsmarkRand : public StableRand
{
public:
    HoltsmarkRand(double scale = 1, double location = 0) : StableRand(1.5, 0.0, scale, location) {}
    std::string Name() const override;
};

/**
 * @brief The LandauRand class
 *
 * X ~ Landau(σ, μ)
 *
 * Related distributions:
 * X ~ Stable(1, 1, σ, μ)
 */
class RANDLIBSHARED_EXPORT LandauRand : public StableRand
{
public:
    LandauRand(double scale = 1, double location = 0) : StableRand(1.0, 1.0, scale, location) {}
    std::string Name() const override;
};

#endif // STABLERAND_H

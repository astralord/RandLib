#ifndef STABLERAND_H
#define STABLERAND_H

#include <functional>
#include "ContinuousDistribution.h"

/**
 * @brief The StableDistribution class <BR>
 * Abstract class for Stable distribution
 *
 * Notation: X ~ S(α, β, γ, μ)
 *
 * Related distributions: <BR>
 * If X ~ Normal(μ, σ), then X ~ S(2, 0, σ/√2, μ) <BR>
 * If X ~ Cauchy(μ, γ), then X ~ S(1, 0, γ, μ) <BR>
 * If +/-X ~ Levy(μ, γ), then X ~ S(0.5, +/-1, γ, μ)
 */
class RANDLIBSHARED_EXPORT StableDistribution : public ContinuousDistribution
{
protected:
    double alpha = 2; ///< characteristic exponent α
    double beta = 0; ///< skewness β
    double mu = 0; ///< location μ
    double gamma = M_SQRT2; ///< scale γ
    double logGamma = 0.5 * M_LN2; ///< log(γ)

private:
    double alphaInv = 0.5; ///< 1/α
    double zeta = 0; ///< ζ = -β * tan(πα/2)
    double omega = 0; ///< ω = log(1 + ζ^2) / (2α)
    double xi = 0; ///< ξ = atan(-ζ) / α;
    double alpha_alpham1 = 2; ///< α / (α - 1)
    double logGammaPi_2 = M_LNPI - 1.5 * M_LN2; ///< log(γπ/2)

    static constexpr double BIG_NUMBER = 1e9; ///< a.k.a. infinity for pdf and cdf calculations
    static constexpr double ALMOST_TWO = 1.99999; ///< parameter used to identify α close to 2

    enum DISTRIBUTION_ID {
        NORMAL, ///< α = 2
        LEVY, ///< α = 0.5, |β| = 1
        CAUCHY, ///< α = 1, β = 0
        UNITY_EXPONENT, ///< α = 1, β ≠ 0
        COMMON ///< the rest
    };

    DISTRIBUTION_ID distributionId = NORMAL; ///< id of distribution (Gaussian by default)

protected:
    double pdfCoef = 0.5 * (M_LN2 + M_LNPI); ///< hashed coefficient for faster pdf calculations
    double pdftailBound = INFINITY; /// boundary k such that for |x| > k we can use pdf tail approximation
    double cdftailBound = INFINITY; /// boundary k such that for |x| > k we can use cdf tail approximation

    StableDistribution(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~StableDistribution() {}

public:
    SUPPORT_TYPE SupportType() const override;
    double MinValue() const override;
    double MaxValue() const override;

protected:
    void SetParameters(double exponent, double skewness, double scale = 1, double location = 0);

public:
    void SetLocation(double location);
    void SetScale(double scale);

    inline double GetExponent() const { return alpha; }
    inline double GetSkewness() const { return beta; }
    inline double GetScale() const { return gamma; }
    inline double GetLocation() const { return mu; }
    inline double GetLogScale() const { return logGamma; }

    /// Probability distribution functions
protected:
    /**
     * @fn pdfNormal
     * @param x
     * @return probability density function of normal distribution
     */
    double pdfNormal(double x) const;
    /**
     * @fn logpdfNormal
     * @param x
     * @return logarithm of probability density function of normal distribution
     */
    double logpdfNormal(double x) const;
    /**
     * @fn pdfCauchy
     * @param x
     * @return probability density function of Cauchy distribution
     */
    double pdfCauchy(double x) const;
    /**
     * @fn logpdfCauchy
     * @param x
     * @return logarithm of probability density function of Cauchy distribution
     */
    double logpdfCauchy(double x) const;
    /**
     * @fn pdfLevy
     * @param x
     * @return probability density function of Levy distribution
     */
    double pdfLevy(double x) const;
    /**
     * @fn logpdfLevy
     * @param x
     * @return logarithm of probability density function of Levy distribution
     */
    double logpdfLevy(double x) const;
private:
    /**
     * @fn fastpdfExponentiation
     * @param u
     * @return exp(u - exp(u)), accelerated by truncation of input u
     */
    static double fastpdfExponentiation(double u);

    /// functions for pdf calculation for α = 1
    double limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const;
    double integrandAuxForUnityExponent(double theta, double xAdj) const;
    double integrandForUnityExponent(double theta, double xAdj) const;
    double pdfForUnityExponent(double x) const;

    /// functions for pdf calculation for α ≠ 1
    DoublePair seriesZeroParams{};
    /**
     * @fn pdfAtZero
     * @return probability density function for x = 0
     */
    double pdfAtZero() const;
    /**
     * @fn pdfSeriesExpansionAtZero
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @param k number of elements in series
     * @return series expansion of probability density function for x near 0
     */
    double pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const;
    /**
     * @fn pdfSeriesExpansionAtInf
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @param k number of elements in series
     * @return series expansion of probability density function for large x
     */
    double pdfSeriesExpansionAtInf(double logX, double xiAdj, int k) const;
    double pdfTaylorExpansionTailNearCauchy(double x) const; /// ~ f(x, α) - f(x, 1)
    double limitCaseForIntegrandAuxForCommonExponent(double theta, double xiAdj) const;
    double integrandAuxForCommonExponent(double theta, double xAdj, double xiAdj) const;
    double integrandForCommonExponent(double theta, double xAdj, double xiAdj) const;
    double pdfForCommonExponent(double x) const;
public:    
    double f(const double & x) const override;
    double logf(const double & x) const override;

    /// Cumulative distribution functions
protected:
    /**
     * @fn cdfNormal
     * @param x
     * @return cumulative distribution function of normal distribution
     */
    double cdfNormal(double x) const;
    /**
     * @fn cdfNormalCompl
     * @param x
     * @return complementary cumulative distribution function of normal distribution
     */
    double cdfNormalCompl(double x) const;
    /**
     * @fn cdfCauchy
     * @param x
     * @return cumulative distribution function of Cauchy distribution
     */
    double cdfCauchy(double x) const;
    /**
     * @fn cdfCauchyCompl
     * @param x
     * @return complementary cumulative distribution function of Cauchy distribution
     */
    double cdfCauchyCompl(double x) const;
    /**
     * @fn cdfLevy
     * @param x
     * @return cumulative distribution function of Levy distribution
     */
    double cdfLevy(double x) const;
    /**
     * @fn cdfLevyCompl
     * @param x
     * @return complementary cumulative distribution function of Levy distribution
     */
    double cdfLevyCompl(double x) const;
private:
    /**
     * @fn fastcdfExponentiation
     * @param u
     * @return exp(-exp(u)), accelerated by truncation of input u
     */
    static double fastcdfExponentiation(double u);
    /**
     * @fn cdfAtZero
     * @param xiAdj adjusted ξ
     * @return cumulative distribution function for x = 0
     */
    double cdfAtZero(double xiAdj) const;
    /**
     * @fn cdfForUnityExponent
     * @param x
     * @return cumulative distribution function for α = 1, β ≠ 0
     */
    double cdfForUnityExponent(double x) const;
    /**
     * @fn cdfSeriesExpansionAtZero
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @param k number of elements in series
     * @return series expansion of cumulative distribution function for x near 0
     */
    double cdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const;
    /**
     * @fn pdfSeriesExpansionAtInf
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @param k number of elements in series
     * @return series expansion of cumulative distribution function for large x
     */
    double cdfSeriesExpansionAtInf(double logX, double xiAdj, int k) const;
    /**
     * @fn cdfIntegralRepresentation
     * @param absXSt absolute value of standardised x
     * @param xiAdj adjusted ξ
     * @return integral representation of cumulative distribution function for common case of α ≠ 1
     */
    double cdfIntegralRepresentation(double logX, double xiAdj) const;
    /**
     * @fn cdfForCommonExponent
     * @param x
     * @return cumulative distribution function for common case of α ≠ 1
     */
    double cdfForCommonExponent(double x) const;
public:
    double F(const double & x) const override;
    double S(const double & x) const override;

    /// Variates
private:
    /**
     * @fn variateForUnityExponent
     * @return variate, generated by algorithm for α = 1, β ≠ 0
     */
    double variateForUnityExponent() const;
    /**
     * @fn variateForCommonExponent
     * @return variate, generated by algorithm for common case of α ≠ 1
     */
    double variateForCommonExponent() const;
    /**
     * @fn variateForExponentEqualOneHalf
     * @return variate, generated by algorithm for special case of α = 0.5
     */
    double variateForExponentEqualOneHalf() const;
public:
    double Variate() const override;
    void Sample(std::vector<double> &outputData) const override;

public:
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
 * @brief The StableRand class <BR>
 * Stable distribution
 */
class RANDLIBSHARED_EXPORT StableRand : public StableDistribution
{
public:
    StableRand(double exponent = 2, double skewness = 0, double scale = 1, double location = 0) : StableDistribution(exponent, skewness, scale, location) {}
    std::string Name() const override;
    using StableDistribution::SetParameters;
};


/**
 * @brief The HoltsmarkRand class <BR>
 * Holtsmark distribution
 *
 * Notation: X ~ Holtsmark(γ, μ)
 *
 * Related distributions:
 * X ~ S(1.5, 0, γ, μ)
 */
class RANDLIBSHARED_EXPORT HoltsmarkRand : public StableDistribution
{
public:
    HoltsmarkRand(double scale = 1, double location = 0) : StableDistribution(1.5, 0.0, scale, location) {}
    std::string Name() const override;
};


/**
 * @brief The LandauRand class <BR>
 * Landau distribution
 *
 * Notation: X ~ Landau(γ, μ)
 *
 * Related distributions:
 * X ~ S(1, 1, γ, μ)
 */
class RANDLIBSHARED_EXPORT LandauRand : public StableDistribution
{
public:
    LandauRand(double scale = 1, double location = 0) : StableDistribution(1.0, 1.0, scale, location) {}
    std::string Name() const override;
};

#endif // STABLERAND_H

#ifndef STABLERAND_H
#define STABLERAND_H

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
template < typename RealType = double>
class RANDLIBSHARED_EXPORT StableDistribution : public ContinuousDistribution<RealType>
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

    enum DISTRIBUTION_TYPE {
        NORMAL, ///< α = 2
        LEVY, ///< α = 0.5, |β| = 1
        CAUCHY, ///< α = 1, β = 0
        UNITY_EXPONENT, ///< α = 1, β ≠ 0
        GENERAL ///< the rest
    };

    DISTRIBUTION_TYPE distributionType = NORMAL; ///< type of distribution (Gaussian by default)

protected:
    double pdfCoef = 0.5 * (M_LN2 + M_LNPI); ///< hashed coefficient for faster pdf calculations
    double pdftailBound = INFINITY; ///< boundary k such that for |x| > k we can use pdf tail approximation
    double cdftailBound = INFINITY; ///< boundary k such that for |x| > k we can use cdf tail approximation

    StableDistribution(double exponent, double skewness, double scale = 1, double location = 0);
    virtual ~StableDistribution() {}

public:
    SUPPORT_TYPE SupportType() const override;
    RealType MinValue() const override;
    RealType MaxValue() const override;

private:
    void parametersVerification(double exponent, double skewness, double scale);
    void setParametersForNormal();
    void setParametersForCauchy();
    void setParametersForLevy();
    void setParametersForUnityExponent();
    void setParametersForGeneralExponent();
    
protected:
    void SetParameters(double exponent, double skewness, double scale = 1, double location = 0);
    
public:
    void SetLocation(double location);
    void SetScale(double scale);

    /**
     * @fn GetExponent
     * @return characteristic exponent α
     */
    inline double GetExponent() const { return alpha; }
    /**
     * @fn GetSkewness
     * @return skewness parameter β
     */
    inline double GetSkewness() const { return beta; }
    /**
     * @fn GetScale
     * @return scale parameter γ
     */
    inline double GetScale() const { return gamma; }
    /**
     * @fn GetLocation
     * @return location parameter μ
     */
    inline double GetLocation() const { return mu; }
    /**
     * @fn GetLogScale
     * @return logarithm of the scale parameter γ
     */
    inline double GetLogScale() const { return logGamma; }

protected:
    /**
     * @fn pdfNormal
     * @param x
     * @return probability density function of normal distribution
     */
    double pdfNormal(RealType x) const;
    /**
     * @fn logpdfNormal
     * @param x
     * @return logarithm of probability density function of normal distribution
     */
    double logpdfNormal(RealType x) const;
    /**
     * @fn pdfCauchy
     * @param x
     * @return probability density function of Cauchy distribution
     */
    double pdfCauchy(RealType x) const;
    /**
     * @fn logpdfCauchy
     * @param x
     * @return logarithm of probability density function of Cauchy distribution
     */
    double logpdfCauchy(RealType x) const;
    /**
     * @fn pdfLevy
     * @param x
     * @return probability density function of Levy distribution
     */
    double pdfLevy(RealType x) const;
    /**
     * @fn logpdfLevy
     * @param x
     * @return logarithm of probability density function of Levy distribution
     */
    double logpdfLevy(RealType x) const;

private:
    /**
     * @fn fastpdfExponentiation
     * @param u
     * @return exp(u - exp(u)), accelerated by truncation of input u
     */
    static double fastpdfExponentiation(double u);

    /**
     * @fn pdfShortTailExpansionForUnityExponent
     * @param logX
     * @return leading term of pdf short tail series expansion for large x, |β| = 1 and α = 1
     */
    double pdfShortTailExpansionForUnityExponent(double x) const;

    /**
     * @fn limitCaseForIntegrandAuxForUnityExponent
     * @param theta
     * @param xAdj
     * @return large values in the case of closeness to extreme points
     */
    double limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const;
    /**
     * @fn integrandAuxForUnityExponent
     * @param theta
     * @param xAdj
     * @return supplementary value for the integrand used in pdf calculations for α = 1
     */
    double integrandAuxForUnityExponent(double theta, double xAdj) const;
    /**
     * @fn integrandForUnityExponent
     * @param theta
     * @param xAdj
     * @return the value of the integrand used for calculations of pdf for α = 1
     */
    double integrandForUnityExponent(double theta, double xAdj) const;
    /**
     * @fn pdfForUnityExponent
     * @param x
     * @return value of probability density function for α = 1
     */
    double pdfForUnityExponent(double x) const;

    DoublePair seriesZeroParams{};

    /**
     * @fn pdfShortTailExpansionForGeneralExponent
     * @param logX
     * @return leading term of pdf short tail series expansion for |β| = 1 and [(large x and α > 1) or (small x and α < 1)]
     */
    double pdfShortTailExpansionForGeneralExponent(double logX) const;
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
     * @return series expansion of probability density function for large x
     */
    double pdfSeriesExpansionAtInf(double logX, double xiAdj) const;
    /**
     * @fn pdfTaylorExpansionTailNearCauchy
     * @param x
     * @return the Taylor approximated difference ~ f(x, α) - f(x, 1)
     */
    double pdfTaylorExpansionTailNearCauchy(double x) const;
    /**
     * @fn limitCaseForIntegrandAuxForGeneralExponent
     * @param theta
     * @param xiAdj
     * @return large values in the case of closeness to extreme points
     */
    double limitCaseForIntegrandAuxForGeneralExponent(double theta, double xiAdj) const;
    /**
     * @fn integrandAuxForGeneralExponent
     * @param theta
     * @param xAdj
     * @param xiAdj
     * @return supplementary value for the integrand used in pdf calculations for α ≠ 1
     */
    double integrandAuxForGeneralExponent(double theta, double xAdj, double xiAdj) const;
    /**
     * @fn integrandFoGeneralExponent
     * @param theta
     * @param xAdj
     * @param xiAdj
     * @return the value of the integrand used for calculations of pdf for α ≠ 1
     */
    double integrandFoGeneralExponent(double theta, double xAdj, double xiAdj) const;
    /**
     * @fn pdfForGeneralExponent
     * @param x
     * @return value of probability density function for α ≠ 1
     */
    double pdfForGeneralExponent(double x) const;
public:    
    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;

protected:
    /**
     * @fn cdfNormal
     * @param x
     * @return cumulative distribution function of normal distribution
     */
    double cdfNormal(RealType x) const;
    /**
     * @fn cdfNormalCompl
     * @param x
     * @return complementary cumulative distribution function of normal distribution
     */
    double cdfNormalCompl(RealType x) const;
    /**
     * @fn cdfCauchy
     * @param x
     * @return cumulative distribution function of Cauchy distribution
     */
    double cdfCauchy(RealType x) const;
    /**
     * @fn cdfCauchyCompl
     * @param x
     * @return complementary cumulative distribution function of Cauchy distribution
     */
    double cdfCauchyCompl(RealType x) const;
    /**
     * @fn cdfLevy
     * @param x
     * @return cumulative distribution function of Levy distribution
     */
    double cdfLevy(RealType x) const;
    /**
     * @fn cdfLevyCompl
     * @param x
     * @return complementary cumulative distribution function of Levy distribution
     */
    double cdfLevyCompl(RealType x) const;
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
     * @return series expansion of cumulative distribution function for large x
     */
    double cdfSeriesExpansionAtInf(double logX, double xiAdj) const;
    /**
     * @fn cdfIntegralRepresentation
     * @param absXSt absolute value of standardised x
     * @param xiAdj adjusted ξ
     * @return integral representation of cumulative distribution function for general case of α ≠ 1
     */
    double cdfIntegralRepresentation(double logX, double xiAdj) const;
    /**
     * @fn cdfForGeneralExponent
     * @param x
     * @return cumulative distribution function for general case of α ≠ 1
     */
    double cdfForGeneralExponent(double x) const;
public:
    double F(const RealType & x) const override;
    double S(const RealType & x) const override;

private:
    /**
     * @fn variateForUnityExponent
     * @return variate, generated by algorithm for α = 1, β ≠ 0
     */
    double variateForUnityExponent() const;
    /**
     * @fn variateForGeneralExponent
     * @return variate, generated by algorithm for general case of α ≠ 1
     */
    double variateForGeneralExponent() const;
    /**
     * @fn variateForExponentEqualOneHalf
     * @return variate, generated by algorithm for special case of α = 0.5
     */
    double variateForExponentEqualOneHalf() const;
public:
    RealType Variate() const override;
    void Sample(std::vector<RealType> &outputData) const override;

public:
    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

protected:
    /**
     * @fn quantileNormal
     * @param p input parameter in the interval (0, 1)
     * @return quantile for Gaussian distribution
     */
    RealType quantileNormal(double p) const;
    /**
     * @fn quantileNormal1m
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for Gaussian distribution
     */
    RealType quantileNormal1m(double p) const;
    /**
     * @fn quantileCauchy
     * @param p input parameter in the interval (0, 1)
     * @return quantile for Cauchy distribution
     */
    RealType quantileCauchy(double p) const;
    /**
     * @fn quantileCauchy1m
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for Cauchy distribution
     */
    RealType quantileCauchy1m(double p) const;
    /**
     * @fn quantileLevy
     * @param p input parameter in the interval (0, 1)
     * @return quantile for Levy distribution
     */
    RealType quantileLevy(double p) const;
    /**
     * @fn quantileLevy1m
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for Levy distribution
     */
    RealType quantileLevy1m(double p) const;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

protected:
    /**
     * @fn cfNormal
     * @param t positive parameter
     * @return characteristic function for Gaussian distribution
     */
    std::complex<double> cfNormal(double t) const;
    /**
     * @fn cfCauchy
     * @param t positive parameter
     * @return characteristic function for Cauchy distribution
     */
    std::complex<double> cfCauchy(double t) const;
    /**
     * @fn cfLevy
     * @param t positive parameter
     * @return characteristic function for Levy distribution
     */
    std::complex<double> cfLevy(double t) const;

private:
    std::complex<double> CFImpl(double t) const override;
};


/**
 * @brief The StableRand class <BR>
 * Stable distribution
 */
template < typename RealType = double>
class RANDLIBSHARED_EXPORT StableRand : public StableDistribution<RealType>
{
public:
    StableRand(double exponent = 2, double skewness = 0, double scale = 1, double location = 0) : StableDistribution<RealType>(exponent, skewness, scale, location) {}
    String Name() const override;
    
    void SetExponent(double location);
    void SetSkewness(double skewness);
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
template < typename RealType = double>
class RANDLIBSHARED_EXPORT HoltsmarkRand : public StableDistribution<RealType>
{
public:
    HoltsmarkRand(double scale = 1, double location = 0) : StableDistribution<RealType>(1.5, 0.0, scale, location) {}
    String Name() const override;
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
template < typename RealType = double>
class RANDLIBSHARED_EXPORT LandauRand : public StableDistribution<RealType>
{
public:
    LandauRand(double scale = 1, double location = 0) : StableDistribution<RealType>(1.0, 1.0, scale, location) {}
    String Name() const override;
};

#endif // STABLERAND_H

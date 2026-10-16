#pragma once

#include <algorithm>
#include <functional>
#include <vector>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * Quantile of the power-law tails of a stable law with α < 2:
 * P(X > x) ~ C (1 + β) γ^α x^(-α) and P(X < -x) ~ C (1 - β) γ^α x^(-α) as x → ∞, with
 * C = Γ(α) sin(πα/2) / π (1/π for α = 1). Geometric stable laws have the same tails.
 * @param p probability of the tail, small
 * @param upperTail true for P(X > x) = p, false for P(X < x) = p
 * @param alpha characteristic exponent
 * @param beta skewness
 * @param gamma scale
 * @param mu location
 * @return asymptotic quantile, or NaN if that tail is not a power law
 */
inline double stableTailQuantile(double p, bool upperTail, double alpha, double beta, double gamma, double mu)
{
    double weight = upperTail ? 1 + beta : 1 - beta;
    if (alpha >= 2 || weight <= 0)
        return NAN;
    double C = (alpha == 1)
                   ? randlib::consts::inv_pi
                   : std::tgamma(alpha) * std::sin(0.5 * randlib::consts::pi * alpha) * randlib::consts::inv_pi;
    double distance = gamma * std::pow(C * weight / p, 1.0 / alpha);
    return upperTail ? mu + distance : mu - distance;
}

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
template <typename RealType>
class StableDistribution : virtual public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

protected:
    double alpha = 2; ///< characteristic exponent α
    double beta = 0; ///< skewness β
    double mu = 0; ///< location μ
    double gamma = randlib::consts::sqrt2; ///< scale γ
    double logGamma = 0.5 * randlib::consts::ln2; ///< log(γ)

private:
    double alphaInv = 0.5; ///< 1/α
    double zeta = 0; ///< ζ = -β * tan(πα/2)
    double omega = 0; ///< ω = log(1 + ζ^2) / (2α)
    double xi = 0; ///< ξ = atan(-ζ) / α;
    double alpha_alpham1 = 2; ///< α / (α - 1)
    double logGammaPi_2 = randlib::consts::lnpi - 1.5 * randlib::consts::ln2; ///< log(π/(2γ)), used for α = 1

    static constexpr double BIG_NUMBER = 1e9; ///< a.k.a. infinity for pdf and cdf calculations
    static constexpr double ALMOST_TWO = 1.99999; ///< parameter used to identify α close to 2

    enum DISTRIBUTION_TYPE
    {
        NORMAL, ///< α = 2
        LEVY, ///< α = 0.5, |β| = 1
        CAUCHY, ///< α = 1, β = 0
        UNITY_EXPONENT, ///< α = 1, β ≠ 0
        GENERAL ///< the rest
    };

    DISTRIBUTION_TYPE distributionType = NORMAL; ///< type of distribution (Gaussian by default)

protected:
    double pdfCoef =
        0.5 * (randlib::consts::ln2 + randlib::consts::lnpi); ///< hashed coefficient for faster pdf calculations
    double pdftailBound = INFINITY; ///< boundary k such that for |x| > k we can use pdf tail approximation
    double cdftailBound = INFINITY; ///< boundary k such that for |x| > k we can use cdf tail approximation

    StableDistribution(double exponent, double skewness, double scale = 1, double location = 0);

public:
    virtual ~StableDistribution() {}

    [[nodiscard]] Support SupportType() const override;
    [[nodiscard]] RealType MinValue() const override;
    [[nodiscard]] RealType MaxValue() const override;

private:
    void parametersVerification(double exponent, double skewness, double scale);
    /**
     * Throws randlib::not_implemented for parameters where the density, the distribution
     * function and everything derived from them are not yet implemented (exponent close
     * to 1 with non-zero skewness); sampling and the characteristic function still work
     */
    void checkDensitySupport() const;
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
     * @return characteristic exponent α
     */
    [[nodiscard]] inline double GetExponent() const noexcept { return alpha; }
    /**
     * @return skewness parameter β
     */
    [[nodiscard]] inline double GetSkewness() const noexcept { return beta; }
    /**
     * @return scale parameter γ
     */
    [[nodiscard]] inline double GetScale() const noexcept { return gamma; }
    /**
     * @return location parameter μ
     */
    [[nodiscard]] inline double GetLocation() const noexcept { return mu; }
    /**
     * @return logarithm of the scale parameter γ
     */
    [[nodiscard]] inline double GetLogScale() const noexcept { return logGamma; }

protected:
    /**
     * @param x
     * @return probability density function of normal distribution
     */
    double pdfNormal(RealType x) const;
    /**
     * @param x
     * @return logarithm of probability density function of normal distribution
     */
    double logpdfNormal(RealType x) const;
    /**
     * @param x
     * @return probability density function of Cauchy distribution
     */
    double pdfCauchy(RealType x) const;
    /**
     * @param x
     * @return logarithm of probability density function of Cauchy distribution
     */
    double logpdfCauchy(RealType x) const;
    /**
     * @param x
     * @return probability density function of Levy distribution
     */
    double pdfLevy(RealType x) const;
    /**
     * @param x
     * @return logarithm of probability density function of Levy distribution
     */
    double logpdfLevy(RealType x) const;

private:
    /**
     * @param u
     * @return exp(u - exp(u)), accelerated by truncation of input u
     */
    static double fastpdfExponentiation(double u);

    /**
     * @param logX
     * @return leading term of pdf short tail series expansion for large x, |β| = 1 and α = 1
     */
    double pdfShortTailExpansionForUnityExponent(double x) const;

    /**
     * @param theta
     * @param xAdj
     * @return large values in the case of closeness to extreme points
     */
    double limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const;
    /**
     * @param theta
     * @param xAdj
     * @return supplementary value for the integrand used in pdf calculations for α = 1
     */
    double integrandAuxForUnityExponent(double theta, double xAdj) const;
    /**
     * @param theta
     * @param xAdj
     * @return the value of the integrand used for calculations of pdf for α = 1
     */
    double integrandForUnityExponent(double theta, double xAdj) const;
    /**
     * @param x
     * @return value of probability density function for α = 1
     */
    double pdfForUnityExponent(double x) const;
    /**
     * @return the point in (a, b) where u, monotone on that interval, changes sign
     */
    static double rootOfMonotone(const std::function<double(double)> &u, double a, double b);
    /**
     * @return integral of the integrand over [a, b], with breakpoints around the root of u
     */
    static double integrateAroundRoot(const std::function<double(double)> &integrand,
                                      const std::function<double(double)> &u, double a, double b, double root);
    /**
     * @param xAdj
     * @return the root of the integrand exponent u(θ) for α = 1
     */
    double peakForUnityExponent(double xAdj) const;

    DoublePair seriesZeroParams{};

    /**
     * @param logX
     * @return leading term of pdf short tail series expansion for |β| = 1 and [(large x and α > 1) or (small x and α < 1)]
     */
    double pdfShortTailExpansionForGeneralExponent(double logX) const;
    /**
     * @param logX
     * @return leading term of the mass between a totally skewed law with α < 1 and its support boundary
     */
    double cdfShortTailExpansion(double logX) const;
    /**
     * @return probability density function for x = 0
     */
    double pdfAtZero() const;
    /**
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @param k number of elements in series
     * @return series expansion of probability density function for x near 0
     */
    double pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const;
    /**
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @return series expansion of probability density function for large x
     */
    double pdfSeriesExpansionAtInf(double logX, double xiAdj) const;
    /**
     * @param x
     * @return the Taylor approximated difference ~ f(x, α) - f(x, 1)
     */
    double pdfTaylorExpansionTailNearCauchy(double x) const;
    /**
     * @param theta
     * @param xiAdj
     * @return large values in the case of closeness to extreme points
     */
    double limitCaseForIntegrandAuxForGeneralExponent(double theta, double xiAdj) const;
    /**
     * @param theta
     * @param xAdj
     * @param xiAdj
     * @return supplementary value for the integrand used in pdf calculations for α ≠ 1
     */
    double integrandAuxForGeneralExponent(double theta, double xAdj, double xiAdj) const;
    /**
     * @param theta
     * @param xAdj
     * @param xiAdj
     * @return the value of the integrand used for calculations of pdf for α ≠ 1
     */
    double integrandFoGeneralExponent(double theta, double xAdj, double xiAdj) const;
    /**
     * @param x
     * @return value of probability density function for α ≠ 1
     */
    double pdfForGeneralExponent(double x) const;

public:
    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;

protected:
    /**
     * @param x
     * @return cumulative distribution function of normal distribution
     */
    double cdfNormal(RealType x) const;
    /**
     * @param x
     * @return complementary cumulative distribution function of normal distribution
     */
    double cdfNormalCompl(RealType x) const;
    /**
     * @param x
     * @return cumulative distribution function of Cauchy distribution
     */
    double cdfCauchy(RealType x) const;
    /**
     * @param x
     * @return complementary cumulative distribution function of Cauchy distribution
     */
    double cdfCauchyCompl(RealType x) const;
    /**
     * @param x
     * @return cumulative distribution function of Levy distribution
     */
    double cdfLevy(RealType x) const;
    /**
     * @param x
     * @return complementary cumulative distribution function of Levy distribution
     */
    double cdfLevyCompl(RealType x) const;

private:
    /**
     * @param u
     * @return exp(-exp(u)), accelerated by truncation of input u
     */
    static double fastcdfExponentiation(double u);
    /**
     * @param u
     * @return 1 - exp(-exp(u)), computed without cancellation
     */
    static double fastccdfExponentiation(double u);
    /**
     * @param xiAdj adjusted ξ
     * @return cumulative distribution function for x = 0
     */
    double cdfAtZero(double xiAdj) const;
    /**
     * @param x
     * @param complement return the survival function instead
     * @return cumulative distribution function (or its complement) for α = 1, β ≠ 0
     */
    double cdfForUnityExponent(double x, bool complement) const;
    /**
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @param k number of elements in series
     * @return series expansion of cumulative distribution function for x near 0
     */
    double cdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const;
    /**
     * @param logX log(x)
     * @param xiAdj adjusted ξ
     * @return series expansion of cumulative distribution function for large x
     */
    double cdfSeriesExpansionAtInf(double logX, double xiAdj) const;
    /**
     * @param absXSt absolute value of standardised x
     * @param xiAdj adjusted ξ
     * @return integral representation of cumulative distribution function for general case of α ≠ 1
     */
    double cdfIntegralRepresentation(double logX, double xiAdj, bool complement) const;
    /**
     * @param x
     * @return cumulative distribution function for general case of α ≠ 1
     */
    double cdfForGeneralExponent(double x) const;
    /**
     * @param x
     * @return survival function for general case of α ≠ 1, computed directly in the upper tail
     */
    double ccdfForGeneralExponent(double x) const;

public:
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

private:
    /**
     * @return variate, generated by algorithm for α = 1, β ≠ 0
     */
    template <std::uniform_random_bit_generator G>
    double variateForUnityExponent(G &g) const;
    /**
     * @return variate, generated by algorithm for general case of α ≠ 1
     */
    template <std::uniform_random_bit_generator G>
    double variateForGeneralExponent(G &g) const;
    /**
     * @return variate, generated by algorithm for special case of α = 0.5
     */
    template <std::uniform_random_bit_generator G>
    double variateForExponentEqualOneHalf(G &g) const;

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

protected:
    /**
     * @param p input parameter in the interval (0, 1)
     * @return quantile for Gaussian distribution
     */
    RealType quantileNormal(double p) const;
    /**
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for Gaussian distribution
     */
    RealType quantileNormal1m(double p) const;
    /**
     * @param p input parameter in the interval (0, 1)
     * @return quantile for Cauchy distribution
     */
    RealType quantileCauchy(double p) const;
    /**
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for Cauchy distribution
     */
    RealType quantileCauchy1m(double p) const;
    /**
     * @param p input parameter in the interval (0, 1)
     * @return quantile for Levy distribution
     */
    RealType quantileLevy(double p) const;
    /**
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for Levy distribution
     */
    RealType quantileLevy1m(double p) const;

private:
    /**
     * Quantile search that starts far tails from their power-law asymptote
     * @param p probability
     * @param upperTail search for S(x) = p instead of F(x) = p
     */
    RealType searchTailQuantile(double p, bool upperTail) const;
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

protected:
    /**
     * @param t positive parameter
     * @return characteristic function for Gaussian distribution
     */
    std::complex<double> cfNormal(double t) const;
    /**
     * @param t positive parameter
     * @return characteristic function for Cauchy distribution
     */
    std::complex<double> cfCauchy(double t) const;
    /**
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
template <typename RealType>
class StableRand final : public StableDistribution<RealType>
{
public:
    StableRand(double exponent = 2, double skewness = 0, double scale = 1, double location = 0)
        : StableDistribution<RealType>(exponent, skewness, scale, location)
    {
    }
    [[nodiscard]] std::string Name() const override;

    void SetExponent(double exponent);
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
template <typename RealType>
class HoltsmarkRand final : public StableDistribution<RealType>
{
public:
    HoltsmarkRand(double scale = 1, double location = 0) : StableDistribution<RealType>(1.5, 0.0, scale, location) {}
    [[nodiscard]] std::string Name() const override;
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
template <typename RealType>
class LandauRand final : public StableDistribution<RealType>
{
public:
    LandauRand(double scale = 1, double location = 0) : StableDistribution<RealType>(1.0, 1.0, scale, location) {}
    [[nodiscard]] std::string Name() const override;
};

} // namespace randlib

#include <randlib/continuous/levy.hpp>
#include <randlib/continuous/cauchy.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/exponential.hpp>

namespace randlib
{

template <typename RealType>
StableDistribution<RealType>::StableDistribution(double exponent, double skewness, double scale, double location)
{
    SetParameters(exponent, skewness, scale, location);
}

template <typename RealType>
Support StableDistribution<RealType>::SupportType() const
{
    if (alpha < 1) {
        if (beta == 1)
            return Support::RightSemiFinite;
        if (beta == -1)
            return Support::LeftSemiFinite;
    }
    return Support::Infinite;
}

template <typename RealType>
RealType StableDistribution<RealType>::MinValue() const
{
    return (alpha < 1 && beta == 1) ? mu : -INFINITY;
}

template <typename RealType>
RealType StableDistribution<RealType>::MaxValue() const
{
    return (alpha < 1 && beta == -1) ? mu : INFINITY;
}

template <typename RealType>
void StableDistribution<RealType>::parametersVerification(double exponent, double skewness, double scale)
{
    if (!(exponent >= 0.1 && exponent <= 2.0))
        throw std::invalid_argument(
            "Stable distribution: exponent should be in the interval [0.1, 2], but it's equal to " +
            std::to_string(exponent));
    if (!(std::fabs(skewness) <= 1.0))
        throw std::invalid_argument(
            "Stable distribution: skewness should be in the interval [-1, 1], but it's equal to " +
            std::to_string(skewness));
    if (!(scale > 0.0))
        throw std::invalid_argument("Stable distribution: scale should be positive, but it's equal to " +
                                    std::to_string(scale));
}

template <typename RealType>
void StableDistribution<RealType>::checkDensitySupport() const
{
    bool exponentNearOne = alpha != 1.0 && std::fabs(alpha - 1.0) < 0.01 && beta != 0.0;
    bool skewnessNearZero = alpha == 1.0 && beta != 0.0 && std::fabs(beta) < 0.01;
    if (exponentNearOne || skewnessNearZero)
        throw randlib::not_implemented(
            "Stable distribution: the density and distribution function are not yet implemented "
            "for exponent within 0.01 of 1 and non-zero skewness, or for exponent 1 and skewness "
            "within 0.01 of 0; sampling and the characteristic function are available");
}

template <typename RealType>
void StableDistribution<RealType>::setParametersForNormal()
{
    distributionType = NORMAL;
    pdfCoef = randlib::consts::ln2 + logGamma + 0.5 * randlib::consts::lnpi;
}

template <typename RealType>
void StableDistribution<RealType>::setParametersForCauchy()
{
    distributionType = CAUCHY;
    pdfCoef = -logGamma - randlib::consts::lnpi;
}

template <typename RealType>
void StableDistribution<RealType>::setParametersForLevy()
{
    distributionType = LEVY;
    pdfCoef = logGamma - randlib::consts::ln2 - randlib::consts::lnpi;
}

template <typename RealType>
void StableDistribution<RealType>::setParametersForUnityExponent()
{
    distributionType = UNITY_EXPONENT;
    pdfCoef = 0.5 / (gamma * std::fabs(beta));
    pdftailBound = 0; // not in the use for now
    /// In S1, X = γZ + (2/π)βγ log γ + μ for α = 1, which is where log γ enters
    logGammaPi_2 = randlib::consts::lnpi - randlib::consts::ln2 - logGamma;
}

template <typename RealType>
void StableDistribution<RealType>::setParametersForGeneralExponent()
{
    distributionType = GENERAL;
    if (beta != 0.0) {
        zeta = -beta * std::tan(randlib::consts::half_pi * alpha);
        omega = 0.5 * alphaInv * std::log1pl(zeta * zeta);
        xi = alphaInv * math::atan(-zeta);
    } else {
        zeta = omega = xi = 0.0;
    }
    pdfCoef = randlib::consts::inv_pi * std::fabs(alpha_alpham1) / gamma;
    pdftailBound = 3.0 / (1.0 + alpha) * randlib::consts::ln10;
    cdftailBound = 3.0 / alpha * randlib::consts::ln10;
    /// define boundaries of region near 0, where we use series expansion
    if (alpha <= ALMOST_TWO) {
        seriesZeroParams.first = std::round(std::min(alpha * alpha * 40 + 1, 10.0));
        seriesZeroParams.second =
            -(alphaInv * 1.5 + 0.5) * randlib::consts::ln10; /// corresponds to interval [10^(-15.5), ~0.056]
    } else {
        seriesZeroParams.first = 85;
        seriesZeroParams.second = randlib::consts::ln2 + randlib::consts::ln3; ///< corresponds to 6
    }
}

template <typename RealType>
void StableDistribution<RealType>::SetParameters(double exponent, double skewness, double scale, double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Stable distribution: location should be finite");
    parametersVerification(exponent, skewness, scale);
    alpha = exponent;
    alphaInv = 1.0 / alpha;
    beta = skewness;
    mu = location;
    gamma = scale;
    logGamma = std::log(gamma);
    alpha_alpham1 = alpha / (alpha - 1.0);

    /// Set id of distribution
    if (alpha == 2.0)
        return setParametersForNormal();
    if (alpha == 1.0)
        return (beta == 0.0) ? setParametersForCauchy() : setParametersForUnityExponent();
    if (alpha == 0.5 && std::fabs(beta) == 1.0)
        return setParametersForLevy();
    return setParametersForGeneralExponent();
}

template <typename RealType>
void StableDistribution<RealType>::SetLocation(double location)
{
    SetParameters(alpha, beta, gamma, location);
}

template <typename RealType>
void StableDistribution<RealType>::SetScale(double scale)
{
    SetParameters(alpha, beta, scale, mu);
}

template <typename RealType>
double StableDistribution<RealType>::pdfNormal(RealType x) const
{
    return std::exp(logpdfNormal(x));
}

template <typename RealType>
double StableDistribution<RealType>::logpdfNormal(RealType x) const
{
    double y = x - mu;
    y *= 0.5 / gamma;
    y *= y;
    y += pdfCoef;
    return -y;
}

template <typename RealType>
double StableDistribution<RealType>::pdfCauchy(RealType x) const
{
    double y = x - mu;
    y *= y;
    y /= gamma;
    y += gamma;
    return randlib::consts::inv_pi / y;
}

template <typename RealType>
double StableDistribution<RealType>::logpdfCauchy(RealType x) const
{
    double x0 = x - mu;
    x0 /= gamma;
    double xSq = x0 * x0;
    return pdfCoef - std::log1pl(xSq);
}

template <typename RealType>
double StableDistribution<RealType>::pdfLevy(RealType x) const
{
    return (x <= mu) ? 0.0 : std::exp(logpdfLevy(x));
}

template <typename RealType>
double StableDistribution<RealType>::logpdfLevy(RealType x) const
{
    double x0 = x - mu;
    if (x0 <= 0.0)
        return -INFINITY;
    double y = gamma / x0;
    y += 3 * std::log(x0);
    y -= pdfCoef;
    return -0.5 * y;
}

template <typename RealType>
double StableDistribution<RealType>::fastpdfExponentiation(double u)
{
    /// exp(u - exp(u)); for very negative u, exp(u) alone, and for large u it underflows to 0 by itself
    return (u < -700) ? std::exp(u) : std::exp(u - std::exp(u));
}

template <typename RealType>
double StableDistribution<RealType>::pdfShortTailExpansionForUnityExponent(double x) const
{
    if (x > 10)
        return 0.0;
    double xm1 = x - 1.0;
    double y = 0.5 * xm1 - std::exp(xm1);
    y -= 0.5 * (randlib::consts::ln2 + randlib::consts::lnpi);
    return std::exp(y);
}

template <typename RealType>
double StableDistribution<RealType>::limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const
{
    if (theta > 0.0) {
        if (beta > 0.0)
            return BIG_NUMBER;
        return (beta == -1) ? xAdj - 1.0 : -BIG_NUMBER;
    }
    if (beta < 0.0)
        return BIG_NUMBER;
    return (beta == 1) ? xAdj - 1.0 : -BIG_NUMBER;
}

template <typename RealType>
double StableDistribution<RealType>::integrandAuxForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi)
        return limitCaseForIntegrandAuxForUnityExponent(theta, xAdj);
    if (theta == 0.0)
        return xAdj + randlib::consts::lnpi - randlib::consts::ln2;
    double thetaAdj = (randlib::consts::half_pi + beta * theta) / std::cos(theta);
    double u = std::log(thetaAdj);
    u += thetaAdj * std::sin(theta) / beta;
    return std::isfinite(u) ? u + xAdj : limitCaseForIntegrandAuxForUnityExponent(theta, xAdj);
}

template <typename RealType>
double StableDistribution<RealType>::integrandForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi)
        return 0.0;
    double u = integrandAuxForUnityExponent(theta, xAdj);
    return fastpdfExponentiation(u);
}

template <typename RealType>
double StableDistribution<RealType>::rootOfMonotone(const std::function<double(double)> &u, double a, double b)
{
    /// Bisection on the sign of u, which changes exactly once on (a, b)
    bool increasing = u(b) > u(a);
    for (int i = 0; i < 200; ++i) {
        double middle = 0.5 * (a + b);
        if (!(middle > a && middle < b))
            break;
        bool belowRoot = (u(middle) < 0) == increasing;
        (belowRoot ? a : b) = middle;
    }
    return 0.5 * (a + b);
}

template <typename RealType>
double StableDistribution<RealType>::integrateAroundRoot(const std::function<double(double)> &integrand,
                                                         const std::function<double(double)> &u, double a, double b,
                                                         double root)
{
    /// In the tails the integrand is concentrated in a region of width ~1/|u'(root)| around the
    /// root of u, which can be far narrower than [a, b]: a quadrature rule over the whole
    /// interval may not place a single node there. Break the interval at log-spaced distances
    /// from the root, scaled by that width.
    double h = 1e-7 * std::max(std::min(root - a, b - root), 1e-300);
    double slope = (u(root + h) - u(root - h)) / (2 * h);
    double width = (std::isfinite(slope) && slope != 0.0) ? 1.0 / std::fabs(slope) : b - a;
    std::vector<double> points = {a, root, b};
    for (double scale : {1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0}) {
        points.push_back(root - scale * width);
        points.push_back(root + scale * width);
    }
    std::sort(points.begin(), points.end());
    double sum = 0.0, left = a;
    for (double point : points) {
        if (point <= left || point > b)
            continue;
        sum += math::integrate(integrand, left, point, 1e-11, 262144);
        left = point;
    }
    return sum;
}

template <typename RealType>
double StableDistribution<RealType>::peakForUnityExponent(double xAdj) const
{
    return rootOfMonotone([this, xAdj](double theta) { return integrandAuxForUnityExponent(theta, xAdj); },
                          -randlib::consts::half_pi, randlib::consts::half_pi);
}

template <typename RealType>
double StableDistribution<RealType>::pdfForUnityExponent(double x) const
{
    double xSt = (x - mu) / gamma;
    double xAdj = -randlib::consts::half_pi * xSt / beta - logGammaPi_2;
    /// The integrand g exp(-g) with g = exp(u) peaks where u = 0, which approaches ±π/2 in the
    /// tails; integrating on both sides of the peak lets the quadrature resolve it
    double theta0 = peakForUnityExponent(xAdj);
    auto integrand = [this, xAdj](double theta) { return integrandForUnityExponent(theta, xAdj); };
    auto u = [this, xAdj](double theta) { return integrandAuxForUnityExponent(theta, xAdj); };
    return pdfCoef * integrateAroundRoot(integrand, u, -randlib::consts::half_pi, randlib::consts::half_pi, theta0);
}

template <typename RealType>
double StableDistribution<RealType>::pdfShortTailExpansionForGeneralExponent(double logX) const
{
    double logAlpha = std::log(alpha);
    double log1mAlpha = (alpha < 1) ? std::log1pl(-alpha) : std::log(alpha - 1);
    double temp = logX - logAlpha;
    double y = std::exp(alpha_alpham1 * temp);
    y *= -std::fabs(1.0 - alpha);
    double z = (1.0 - 0.5 * alpha) / (alpha - 1) * temp;
    z -= 0.5 * (randlib::consts::ln2 + randlib::consts::lnpi + logAlpha + log1mAlpha);
    return std::exp(y + z);
}

template <typename RealType>
double StableDistribution<RealType>::cdfShortTailExpansion(double logX) const
{
    /// Mass between a totally skewed law (α < 1) and its support boundary at distance x:
    /// exp(-(1 - α) u) / sqrt(2π α (1 - α) u) with u = (x/α)^(-α/(1-α)), to leading order
    double logU = alpha / (alpha - 1.0) * (logX - std::log(alpha));
    double u = std::exp(logU);
    double y = -(1.0 - alpha) * u -
               0.5 * (randlib::consts::ln2 + randlib::consts::lnpi + std::log(alpha) + std::log1p(-alpha) + logU);
    return std::exp(y);
}

template <typename RealType>
double StableDistribution<RealType>::pdfAtZero() const
{
    double y0 = 0.0;
    if (beta == 0.0)
        y0 = std::tgammal(alphaInv);
    else {
        y0 = std::lgammal(alphaInv) - omega;
        y0 = std::exp(y0) * std::cos(xi);
    }
    return y0 * randlib::consts::inv_pi / alpha;
}

template <typename RealType>
double StableDistribution<RealType>::pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
{
    /// Calculate first term of the sum
    /// (if x is 0, only this term is non-zero)
    double y0 = pdfAtZero();
    double sum = 0.0;
    if (beta == 0.0) {
        /// Symmetric distribution
        for (int n = 1; n <= k; ++n) {
            int n2 = n + n;
            double term = std::lgammal((n2 + 1) / alpha);
            term += n2 * logX;
            term -= math::lfact(n2);
            term = std::exp(term);
            sum += (n & 1) ? -term : term;
        }
    } else {
        /// Asymmetric distribution
        double rhoPi_alpha = randlib::consts::half_pi + xiAdj;
        for (int n = 1; n <= k; ++n) {
            int np1 = n + 1;
            double term = std::lgammal(np1 * alphaInv);
            term += n * logX;
            term -= math::lfact(n);
            term = std::exp(term - omega);
            term *= std::sin(np1 * rhoPi_alpha);
            sum += (n & 1) ? -term : term;
        }
    }
    return y0 + sum * randlib::consts::inv_pi / alpha;
}

template <typename RealType>
double StableDistribution<RealType>::pdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; ///< number of elements in the series
    double rhoPi = randlib::consts::half_pi + xiAdj;
    rhoPi *= alpha;
    double sum = 0.0;
    for (int n = 1; n <= k; ++n) {
        double aux = n * alpha + 1.0;
        double term = std::lgammal(aux);
        term -= aux * logX;
        term -= math::lfact(n);
        term = std::exp(term - omega);
        term *= std::sin(rhoPi * n);
        sum += (n & 1) ? term : -term;
    }
    return randlib::consts::inv_pi * sum;
}

template <typename RealType>
double StableDistribution<RealType>::pdfTaylorExpansionTailNearCauchy(double x) const
{
    double xSq = x * x;
    double y = 1.0 + xSq;
    double ySq = y * y;
    double z = math::atan(x);
    double zSq = z * z;
    double logY = std::log1pl(xSq);
    double alpham1 = alpha - 1.0;
    double temp = 1.0 - randlib::consts::egamma - 0.5 * logY;
    /// first derivative
    double f_a = temp;
    f_a *= xSq - 1.0;
    f_a += 2 * x * z;
    f_a /= ySq;
    static constexpr long double M_PI_SQ_6 = 1.64493406684822643647l; /// π^2 / 6
    /// second derivative
    double f_aa1 = M_PI_SQ_6;
    f_aa1 += temp * temp;
    f_aa1 -= 1.0 + z * z;
    f_aa1 *= xSq * xSq - 6.0 * xSq + 1.0;
    double f_aa2 = 0.5 + temp;
    f_aa2 *= z;
    f_aa2 *= 8 * x * (xSq - 1.0);
    double f_aa3 = (1.0 - 3 * xSq) * temp;
    f_aa3 -= x * y * z;
    f_aa3 += f_aa3;
    double f_aa = f_aa1 + f_aa2 + f_aa3;
    f_aa /= std::pow(y, 3);
    /// Hashed values of special functions for x = 2, 3, 4
    /// Gamma(x)
    static constexpr int gammaTable[] = {1, 2, 6};
    /// Gamma'(x)
    static constexpr long double gammaDerTable[] = {1.0 - randlib::consts::egamma, 3.0 - 2.0 * randlib::consts::egamma,
                                                    11.0 - 6.0 * randlib::consts::egamma};
    /// Gamma''(x)
    static constexpr long double gammaSecDerTable[] = {0.82368066085287938958l, 2.49292999190269305794l,
                                                       11.1699273161019477314l};
    /// Digamma(x)
    static constexpr long double digammaTable[] = {1.0 - randlib::consts::egamma, 1.5 - randlib::consts::egamma,
                                                   11.0 / 6 - randlib::consts::egamma};
    /// Digamma'(x)
    static constexpr long double digammaDerTable[] = {M_PI_SQ_6 - 1.0, M_PI_SQ_6 - 1.25, M_PI_SQ_6 - 49.0 / 36};
    /// Digamma''(x)
    static constexpr long double digammaSecDerTable[] = {-0.40411380631918857080l, -0.15411380631918857080l,
                                                         -0.08003973224511449673l};
    /// third derivative
    double gTable[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        double g_11 = 0.25 * gammaTable[i] * logY * logY;
        g_11 -= gammaDerTable[i] * logY;
        g_11 += gammaSecDerTable[i];
        double aux = digammaTable[i] - 0.5 * logY;
        double g_12 = aux;
        double zip2 = z * (i + 2);
        double cosZNu = std::cos(zip2), zSinZNu = z * std::sin(zip2);
        g_12 *= cosZNu;
        g_12 -= zSinZNu;
        double g_1 = g_11 * g_12;
        double g_21 = -gammaTable[i] * logY + 2 * gammaDerTable[i];
        double g_22 = -zSinZNu * aux;
        g_22 -= zSq * cosZNu;
        g_22 += cosZNu * digammaDerTable[i];
        double g_2 = g_21 * g_22;
        double g_3 = -zSq * cosZNu * aux;
        g_3 -= 2 * zSinZNu * digammaDerTable[i];
        g_3 += zSq * zSinZNu;
        g_3 += cosZNu * digammaSecDerTable[i];
        g_3 *= gammaTable[i];
        double g = g_1 + g_2 + g_3;
        g *= std::pow(y, -0.5 * i - 1);
        gTable[i] = g;
    }
    double f_aaa = -gTable[0] + 3 * gTable[1] - gTable[2];
    /// sum all three derivatives
    double tail = f_a * alpham1;
    tail += 0.5 * f_aa * alpham1 * alpham1;
    tail += std::pow(alpham1, 3) * f_aaa / 6.0;
    tail /= randlib::consts::pi;
    return tail;
}

template <typename RealType>
double StableDistribution<RealType>::limitCaseForIntegrandAuxForGeneralExponent(double theta, double xiAdj) const
{
    /// We got numerical error, need to investigate to which extreme point we are closer
    if (theta < 0.5 * (randlib::consts::half_pi - xiAdj))
        return alpha < 1 ? -BIG_NUMBER : BIG_NUMBER;
    return alpha < 1 ? BIG_NUMBER : -BIG_NUMBER;
}

template <typename RealType>
double StableDistribution<RealType>::integrandAuxForGeneralExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi || theta <= -xiAdj)
        return limitCaseForIntegrandAuxForGeneralExponent(theta, xiAdj);
    double thetaAdj = alpha * (theta + xiAdj);
    double sinThetaAdj = std::sin(thetaAdj);
    double y = std::log(std::cos(theta));
    y -= alpha * std::log(sinThetaAdj);
    y /= alpha - 1.0;
    y += std::log(std::cos(thetaAdj - theta));
    y += xAdj;
    return std::isfinite(y) ? y : limitCaseForIntegrandAuxForGeneralExponent(theta, xiAdj);
}

template <typename RealType>
double StableDistribution<RealType>::integrandFoGeneralExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi)
        return 0.0;
    if (theta <= -xiAdj)
        return 0.0;
    double u = integrandAuxForGeneralExponent(theta, xAdj, xiAdj);
    return fastpdfExponentiation(u);
}

template <typename RealType>
double StableDistribution<RealType>::pdfForGeneralExponent(double x) const
{
    /// Standardize
    double xSt = (x - mu) / gamma;
    double absXSt = xSt;
    /// +- xi
    double xiAdj = xi;
    if (xSt > 0) {
        if (alpha < 1 && beta == -1)
            return 0.0;
    } else {
        if (alpha < 1 && beta == 1)
            return 0.0;
        absXSt = -xSt;
        xiAdj = -xi;
    }

    /// If α is too close to 1 and distribution is symmetric, then we approximate using Taylor series
    if (beta == 0.0 && std::fabs(alpha - 1.0) < 0.01)
        return pdfCauchy(x) + pdfTaylorExpansionTailNearCauchy(absXSt) / gamma;

    /// If x = 0, we know the analytic solution
    if (xSt == 0.0)
        return pdfAtZero() / gamma;

    double logAbsX = std::log(absXSt) - omega;

    /// If x is too close to 0, we do series expansion avoiding numerical problems
    if (logAbsX < seriesZeroParams.second) {
        if (alpha < 1 && std::fabs(beta) == 1)
            return pdfShortTailExpansionForGeneralExponent(logAbsX);
        return pdfSeriesExpansionAtZero(logAbsX, xiAdj, seriesZeroParams.first) / gamma;
    }

    /// If x is large enough we use tail approximation
    if (logAbsX > pdftailBound && alpha <= ALMOST_TWO) {
        /// for α > 1 and |β| = 1 only the tail on the side of -β is short; the other one is heavy
        double skewnessTowardsX = (xSt > 0) ? beta : -beta;
        if (alpha > 1 && skewnessTowardsX == -1)
            return pdfShortTailExpansionForGeneralExponent(logAbsX);
        return pdfSeriesExpansionAtInf(logAbsX, xiAdj) / gamma;
    }

    double xAdj = alpha_alpham1 * logAbsX;

    /// Integrate around the peak of the integrand, where its exponent u vanishes
    auto u = [this, xAdj, xiAdj](double theta) { return integrandAuxForGeneralExponent(theta, xAdj, xiAdj); };
    auto integrand = [this, xAdj, xiAdj](double theta) { return integrandFoGeneralExponent(theta, xAdj, xiAdj); };
    double theta0 = rootOfMonotone(u, -xiAdj, randlib::consts::half_pi);
    double integral = integrateAroundRoot(integrand, u, -xiAdj, randlib::consts::half_pi, theta0);
    double res = pdfCoef * integral / absXSt;

    /// Finally we check if α is not too close to 2
    if (alpha <= ALMOST_TWO)
        return res;

    /// If α is near 2, we use tail aprroximation for large x
    /// and compare it with integral representation
    double alphap1 = alpha + 1.0;
    double tail = std::lgammal(alphap1);
    tail -= alphap1 * logAbsX;
    tail = std::exp(tail);
    tail *= (1.0 - 0.5 * alpha) / gamma;
    return std::max(tail, res);
}

template <typename RealType>
double StableDistribution<RealType>::f(const RealType &x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return pdfNormal(x);
    case CAUCHY:
        return pdfCauchy(x);
    case LEVY:
        return (beta > 0) ? pdfLevy(x) : pdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return pdfForUnityExponent(x);
    case GENERAL:
        return pdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template <typename RealType>
double StableDistribution<RealType>::logf(const RealType &x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return logpdfNormal(x);
    case CAUCHY:
        return logpdfCauchy(x);
    case LEVY:
        return (beta > 0) ? logpdfLevy(x) : logpdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return std::log(pdfForUnityExponent(x));
    case GENERAL:
        return std::log(pdfForGeneralExponent(x));
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template <typename RealType>
double StableDistribution<RealType>::cdfNormal(RealType x) const
{
    double y = mu - x;
    y *= 0.5 / gamma;
    return 0.5 * std::erfc(y);
}

template <typename RealType>
double StableDistribution<RealType>::cdfNormalCompl(RealType x) const
{
    double y = x - mu;
    y *= 0.5 / gamma;
    return 0.5 * std::erfc(y);
}

template <typename RealType>
double StableDistribution<RealType>::cdfCauchy(RealType x) const
{
    double x0 = x - mu;
    x0 /= gamma;
    /// In the lower tail 1/2 + atan(x)/π cancels; atan(-1/x)/π does not
    return (x0 < -1.0) ? randlib::consts::inv_pi * std::atan(-1.0 / x0) : 0.5 + randlib::consts::inv_pi * std::atan(x0);
}

template <typename RealType>
double StableDistribution<RealType>::cdfCauchyCompl(RealType x) const
{
    double x0 = mu - x;
    x0 /= gamma;
    return (x0 < -1.0) ? randlib::consts::inv_pi * std::atan(-1.0 / x0) : 0.5 + randlib::consts::inv_pi * std::atan(x0);
}

template <typename RealType>
double StableDistribution<RealType>::cdfLevy(RealType x) const
{
    if (x <= mu)
        return 0;
    double y = x - mu;
    y += y;
    y = gamma / y;
    y = std::sqrt(y);
    return std::erfc(y);
}

template <typename RealType>
double StableDistribution<RealType>::cdfLevyCompl(RealType x) const
{
    if (x <= mu)
        return 1.0;
    double y = x - mu;
    y += y;
    y = gamma / y;
    y = std::sqrt(y);
    return std::erf(y);
}

template <typename RealType>
double StableDistribution<RealType>::fastcdfExponentiation(double u)
{
    /// exp(-exp(u)), which underflows to 0 by itself for large u
    return std::exp(-std::exp(u));
}

template <typename RealType>
double StableDistribution<RealType>::fastccdfExponentiation(double u)
{
    return (u > 5.0) ? 1.0 : -std::expm1(-std::exp(u));
}

template <typename RealType>
double StableDistribution<RealType>::cdfAtZero(double xiAdj) const
{
    return 0.5 - randlib::consts::inv_pi * xiAdj;
}

template <typename RealType>
double StableDistribution<RealType>::cdfForUnityExponent(double x, bool complement) const
{
    /// y = (1/π) ∫ exp(-exp(u(θ))) dθ is F for β > 0 and S for β < 0; the complementary
    /// integrand 1 - exp(-exp(u)) gives the other one without cancellation. The integrand
    /// switches between 0 and 1 around the root θ0 of u, so integrate on both sides of it.
    double xSt = (x - mu) / gamma;
    double xAdj = -randlib::consts::half_pi * xSt / beta - logGammaPi_2;
    bool useComplement = complement == (beta > 0);
    auto integrand = [this, xAdj, useComplement](double theta) {
        double u = integrandAuxForUnityExponent(theta, xAdj);
        return useComplement ? fastccdfExponentiation(u) : fastcdfExponentiation(u);
    };
    double theta0 = peakForUnityExponent(xAdj);
    auto u = [this, xAdj](double theta) { return integrandAuxForUnityExponent(theta, xAdj); };
    return randlib::consts::inv_pi *
           integrateAroundRoot(integrand, u, -randlib::consts::half_pi, randlib::consts::half_pi, theta0);
}

template <typename RealType>
double StableDistribution<RealType>::cdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
{
    /// Calculate first term of the sum
    /// (if x = 0, only this term is non-zero)
    double y0 = cdfAtZero(xiAdj);
    double sum = 0.0;
    if (beta == 0.0) {
        /// Symmetric distribution
        for (int m = 0; m <= k; ++m) {
            int m2p1 = 2 * m + 1;
            double term = std::lgammal(m2p1 * alphaInv);
            term += m2p1 * logX;
            term -= math::lfact(m2p1);
            term = std::exp(term);
            sum += (m & 1) ? -term : term;
        }
    } else {
        /// Asymmetric distribution
        double rhoPi_alpha = randlib::consts::half_pi + xiAdj;
        for (int n = 1; n <= k; ++n) {
            double term = std::lgammal(n * alphaInv);
            term += n * logX;
            term -= math::lfact(n);
            term = std::exp(term);
            term *= std::sin(n * rhoPi_alpha);
            sum += (n & 1) ? term : -term;
        }
    }
    return y0 + sum * randlib::consts::inv_pi * alphaInv;
}

template <typename RealType>
double StableDistribution<RealType>::cdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; /// number of elements in the series
    double rhoPi = randlib::consts::half_pi + xiAdj;
    rhoPi *= alpha;
    double sum = 0.0;
    for (int n = 1; n <= k; ++n) {
        double aux = n * alpha;
        double term = std::lgammal(aux);
        term -= aux * logX;
        term -= math::lfact(n);
        term = std::exp(term);
        term *= std::sin(rhoPi * n);
        sum += (n & 1) ? term : -term;
    }
    return randlib::consts::inv_pi * sum;
}

template <typename RealType>
double StableDistribution<RealType>::cdfIntegralRepresentation(double logX, double xiAdj, bool complement) const
{
    /// (1/π) ∫_{-ξ}^{π/2} exp(-exp(u(θ))) dθ, or the same integral of 1 - exp(-exp(u(θ)))
    double xAdj = alpha_alpham1 * logX;
    auto u = [this, xAdj, xiAdj](double theta) { return integrandAuxForGeneralExponent(theta, xAdj, xiAdj); };
    auto integrand = [&u, complement](double theta) {
        return complement ? fastccdfExponentiation(u(theta)) : fastcdfExponentiation(u(theta));
    };
    double theta0 = rootOfMonotone(u, -xiAdj, randlib::consts::half_pi);
    return randlib::consts::inv_pi * integrateAroundRoot(integrand, u, -xiAdj, randlib::consts::half_pi, theta0);
}

template <typename RealType>
double StableDistribution<RealType>::cdfForGeneralExponent(double x) const
{
    double xSt = (x - mu) / gamma; /// Standardize
    if (xSt == 0) {
        /// totally skewed laws with α < 1 start (or end) exactly at μ
        if (alpha < 1 && std::fabs(beta) == 1)
            return (beta > 0) ? 0.0 : 1.0;
        return cdfAtZero(xi);
    }
    if (xSt > 0.0) {
        if (alpha > 1.0)
            return 1.0 - ccdfForGeneralExponent(x);
        /// α < 1: F is small near 0, so compute it directly rather than as 1 - S
        if (beta == -1.0)
            return 1.0;
        double logAbsX = std::log(xSt) - omega;
        if (logAbsX > cdftailBound)
            return 1.0 - cdfSeriesExpansionAtInf(logAbsX, xi);
        if (logAbsX < seriesZeroParams.second)
            return (beta == 1.0) ? cdfShortTailExpansion(logAbsX)
                                 : cdfSeriesExpansionAtZero(logAbsX, xi, seriesZeroParams.first);
        double atZero = (beta == 1.0) ? 0.0 : cdfAtZero(xi);
        return atZero + cdfIntegralRepresentation(logAbsX, xi, false);
    }
    /// For x < 0 we use relation F(x, ξ) = S(-x, -ξ)
    double logAbsX = std::log(-xSt) - omega;
    if (logAbsX < seriesZeroParams.second) {
        if (alpha < 1.0 && beta == -1.0)
            return 1.0 - cdfShortTailExpansion(logAbsX);
        return 1.0 - cdfSeriesExpansionAtZero(logAbsX, -xi, seriesZeroParams.first);
    }
    if (logAbsX > cdftailBound)
        return cdfSeriesExpansionAtInf(logAbsX, -xi);
    if (alpha > 1.0)
        return cdfIntegralRepresentation(logAbsX, -xi, false);
    return (beta == 1.0) ? 0.0 : cdfIntegralRepresentation(logAbsX, -xi, true);
}

template <typename RealType>
double StableDistribution<RealType>::ccdfForGeneralExponent(double x) const
{
    double xSt = (x - mu) / gamma; /// Standardize
    if (xSt <= 0)
        return 1.0 - cdfForGeneralExponent(x);
    double logAbsX = std::log(xSt) - omega;
    /// If x is too close to 0, we do series expansion avoiding numerical problems
    if (logAbsX < seriesZeroParams.second) {
        if (alpha < 1.0 && beta == 1.0)
            return 1.0 - cdfShortTailExpansion(logAbsX);
        if (alpha < 1.0 && beta == -1.0)
            return cdfShortTailExpansion(logAbsX);
        return 1.0 - cdfSeriesExpansionAtZero(logAbsX, xi, seriesZeroParams.first);
    }
    /// If x is large enough we use tail approximation
    if (logAbsX > cdftailBound)
        return cdfSeriesExpansionAtInf(logAbsX, xi);
    if (alpha > 1.0)
        return cdfIntegralRepresentation(logAbsX, xi, false);
    return (beta == -1.0) ? 0.0 : cdfIntegralRepresentation(logAbsX, xi, true);
}

template <typename RealType>
double StableDistribution<RealType>::F(const RealType &x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return cdfNormal(x);
    case CAUCHY:
        return cdfCauchy(x);
    case LEVY:
        return (beta > 0) ? cdfLevy(x) : cdfLevyCompl(2 * mu - x);
    case UNITY_EXPONENT:
        return cdfForUnityExponent(x, false);
    case GENERAL:
        return cdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template <typename RealType>
double StableDistribution<RealType>::S(const RealType &x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return cdfNormalCompl(x);
    case CAUCHY:
        return cdfCauchyCompl(x);
    case LEVY:
        return (beta > 0) ? cdfLevyCompl(x) : cdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return cdfForUnityExponent(x, true);
    case GENERAL:
        return ccdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
double StableDistribution<RealType>::variateForUnityExponent(G &g) const
{
    RealType U = randlib::consts::pi * detail::uniform01<RealType>(g) - randlib::consts::half_pi;
    RealType W = detail::standardExponential<RealType>(g);
    RealType pi_2pBetaU = randlib::consts::half_pi + beta * U;
    RealType Y = W * std::cos(U) / pi_2pBetaU;
    RealType X = std::log(Y);
    X += logGammaPi_2;
    X *= -beta;
    X += pi_2pBetaU * std::tan(U);
    X *= randlib::consts::two_inv_pi;
    return mu + gamma * X;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
double StableDistribution<RealType>::variateForGeneralExponent(G &g) const
{
    RealType U = randlib::consts::pi * detail::uniform01<RealType>(g) - randlib::consts::half_pi;
    RealType W = detail::standardExponential<RealType>(g);
    RealType alphaUpxi = alpha * (U + xi);
    RealType X = std::sin(alphaUpxi);
    RealType W_adj = W / std::cos(U - alphaUpxi);
    X *= W_adj;
    RealType R = omega - alphaInv * std::log(W_adj * std::cos(U));
    X *= std::exp(R);
    return mu + gamma * X;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
double StableDistribution<RealType>::variateForExponentEqualOneHalf(G &g) const
{
    RealType Z1 = detail::standardNormal<RealType>(g);
    RealType Z2 = detail::standardNormal<RealType>(g);
    RealType temp1 = (1.0 + beta) / Z1, temp2 = (1.0 - beta) / Z2;
    RealType var = temp1 - temp2;
    var *= temp1 + temp2;
    var *= 0.25;
    return mu + gamma * var;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType StableDistribution<RealType>::operator()(G &g) const
{
    switch (distributionType) {
    case NORMAL:
        return mu + randlib::consts::sqrt2 * gamma * detail::standardNormal<RealType>(g);
    case CAUCHY:
        return mu + gamma * CauchyRand<RealType>::StandardVariate(g);
    case LEVY:
        return mu + math::sign(beta) * gamma * LevyRand<RealType>::StandardVariate(g);
    case UNITY_EXPONENT:
        return variateForUnityExponent(g);
    case GENERAL:
        return (alpha == 0.5) ? variateForExponentEqualOneHalf(g) : variateForGeneralExponent(g);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template <typename RealType>
void StableDistribution<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    switch (distributionType) {
    case NORMAL: {
        double stdev = randlib::consts::sqrt2 * gamma;
        for (RealType &var : outputData)
            var = mu + stdev * detail::standardNormal<RealType>(g);
    } break;
    case CAUCHY: {
        for (RealType &var : outputData)
            var = mu + gamma * CauchyRand<RealType>::StandardVariate(g);
    } break;
    case LEVY: {
        if (beta > 0) {
            for (RealType &var : outputData)
                var = mu + gamma * LevyRand<RealType>::StandardVariate(g);
        } else {
            for (RealType &var : outputData)
                var = mu - gamma * LevyRand<RealType>::StandardVariate(g);
        }
    } break;
    case UNITY_EXPONENT: {
        for (RealType &var : outputData)
            var = variateForUnityExponent(g);
    } break;
    case GENERAL: {
        if (alpha == 0.5) {
            for (RealType &var : outputData)
                var = variateForExponentEqualOneHalf(g);
        } else {
            for (RealType &var : outputData)
                var = variateForGeneralExponent(g);
        }
    } break;
    default:
        break;
    }
}

template <typename RealType>
long double StableDistribution<RealType>::Mean() const
{
    if (alpha > 1)
        return mu;
    if (beta == 1)
        return INFINITY;
    return (beta == -1) ? -INFINITY : NAN;
}

template <typename RealType>
long double StableDistribution<RealType>::Variance() const
{
    return (distributionType == NORMAL) ? 2 * gamma * gamma : INFINITY;
}

template <typename RealType>
RealType StableDistribution<RealType>::Mode() const
{
    checkDensitySupport();
    /// For symmetric and normal distributions mode is μ
    if (beta == 0 || distributionType == NORMAL)
        return mu;
    if (distributionType == LEVY)
        return mu + beta * gamma / 3.0;
    return ContinuousDistribution<RealType>::Mode();
}

template <typename RealType>
RealType StableDistribution<RealType>::Median() const
{
    checkDensitySupport();
    /// For symmetric and normal distributions mode is μ
    if (beta == 0 || distributionType == NORMAL)
        return mu;
    return ContinuousDistribution<RealType>::Median();
}

template <typename RealType>
long double StableDistribution<RealType>::Skewness() const
{
    return (distributionType == NORMAL) ? 0 : NAN;
}

template <typename RealType>
long double StableDistribution<RealType>::ExcessKurtosis() const
{
    return (distributionType == NORMAL) ? 0 : NAN;
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileNormal(double p) const
{
    return mu - 2 * gamma * math::erfcinv(2 * p);
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileNormal1m(double p) const
{
    return mu + 2 * gamma * math::erfcinv(2 * p);
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileCauchy(double p) const
{
    return mu - gamma / std::tan(randlib::consts::pi * p);
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileCauchy1m(double p) const
{
    return mu + gamma / std::tan(randlib::consts::pi * p);
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileLevy(double p) const
{
    double y = math::erfcinv(p);
    return mu + 0.5 * gamma / (y * y);
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileLevy1m(double p) const
{
    double y = math::erfinv(p);
    return mu + 0.5 * gamma / (y * y);
}

template <typename RealType>
RealType StableDistribution<RealType>::searchTailQuantile(double p, bool upperTail) const
{
    double tail = std::min(p, 1.0 - p);
    double guess = (tail < 0.05) ? stableTailQuantile(tail, upperTail != (p > 0.5), alpha, beta, gamma, mu) : NAN;
    return this->searchQuantile(p, upperTail, guess, 0.05);
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileImpl(double p) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return quantileNormal(p);
    case CAUCHY:
        return quantileCauchy(p);
    case LEVY:
        return (beta > 0) ? quantileLevy(p) : 2 * mu - quantileLevy1m(p);
    default:
        return searchTailQuantile(p, false);
    }
}

template <typename RealType>
RealType StableDistribution<RealType>::quantileImpl1m(double p) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return quantileNormal1m(p);
    case CAUCHY:
        return quantileCauchy1m(p);
    case LEVY:
        return (beta > 0) ? quantileLevy1m(p) : 2 * mu - quantileLevy(p);
    default:
        return searchTailQuantile(p, true);
    }
}

template <typename RealType>
std::complex<double> StableDistribution<RealType>::cfNormal(double t) const
{
    double gammaT = gamma * t;
    std::complex<double> y(-gammaT * gammaT, mu * t);
    return std::exp(y);
}

template <typename RealType>
std::complex<double> StableDistribution<RealType>::cfCauchy(double t) const
{
    std::complex<double> y(-gamma * t, mu * t);
    return std::exp(y);
}

template <typename RealType>
std::complex<double> StableDistribution<RealType>::cfLevy(double t) const
{
    std::complex<double> y(0.0, -2 * gamma * t);
    y = -std::sqrt(y);
    y += std::complex<double>(0.0, mu * t);
    return std::exp(y);
}

template <typename RealType>
std::complex<double> StableDistribution<RealType>::CFImpl(double t) const
{
    double x = 0;
    switch (distributionType) {
    case NORMAL:
        return cfNormal(t);
    case CAUCHY:
        return cfCauchy(t);
    case LEVY: {
        std::complex<double> phi = cfLevy(t);
        return (beta > 0) ? phi : std::conj(phi);
    }
    case UNITY_EXPONENT:
        x = beta * randlib::consts::two_inv_pi * std::log(t);
        break;
    default:
        x = zeta;
    }
    double re = std::pow(gamma * t, alpha);
    std::complex<double> psi = std::complex<double>(re, re * x - mu * t);
    return std::exp(-psi);
}
template <typename RealType>
std::string StableRand<RealType>::Name() const
{
    return "Stable(" + this->toStringWithPrecision(this->GetExponent()) + ", " +
           this->toStringWithPrecision(this->GetSkewness()) + ", " + this->toStringWithPrecision(this->GetScale()) +
           ", " + this->toStringWithPrecision(this->GetLocation()) + ")";
}

template <typename RealType>
void StableRand<RealType>::SetExponent(double exponent)
{
    this->SetParameters(exponent, this->GetSkewness(), this->GetScale(), this->GetLocation());
}

template <typename RealType>
void StableRand<RealType>::SetSkewness(double skewness)
{
    this->SetParameters(this->GetExponent(), skewness, this->GetScale(), this->GetLocation());
}
template <typename RealType>
std::string HoltsmarkRand<RealType>::Name() const
{
    return "Holtsmark(" + this->toStringWithPrecision(this->GetScale()) + ", " +
           this->toStringWithPrecision(this->GetLocation()) + ")";
}
template <typename RealType>
std::string LandauRand<RealType>::Name() const
{
    return "Landau(" + this->toStringWithPrecision(this->GetScale()) + ", " +
           this->toStringWithPrecision(this->GetLocation()) + ")";
}

} // namespace randlib

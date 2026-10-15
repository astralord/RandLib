#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/stable.hpp>

namespace randlib
{

/**
 * @brief The GeneralGeometricStableDistribution class
 * Abstract class for general interface of Asymmetric Laplace
 * and Geometric-Stable distributions
 */
template <typename RealType>
class GeneralGeometricStableDistribution : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    StableRand<RealType> Z{};

protected:
    double alpha = 2; ///< characteristic exponent α
    double alphaInv = 0.5; /// 1/α
    double beta = 0; ///< skewness β
    double mu = 0; ///< location μ
    double m = 0; ///< shift m
    double gamma = randlib::consts::sqrt2; ///< scale γ
    double logGamma = 0.5 * randlib::consts::ln2; ///< log(γ)

    /// parameters for α = 2
    double kappa = 1; ///< asymmetry coefficient κ
    double kappaInv = 1; ///< 1 / κ
    double kappaSq = 1; ///< κ^2
    double log1pKappaSq = randlib::consts::ln2; ///< log(1 + κ^2)
    double logKappa = 0; /// log(κ)

    GeneralGeometricStableDistribution(double exponent, double skewness, double scale = 1.0, double location = 0.0,
                                       double shift = 0.0);

public:
    virtual ~GeneralGeometricStableDistribution() {}

protected:
    void SetParameters(double exponent, double skewness, double scale = 1.0, double location = 0.0, double shift = 0.);
    void SetLocation(double location);
    void SetShift(double shift);
    void SetScale(double scale);
    void SetAsymmetry(double asymmetry);

    enum DISTRIBUTION_TYPE
    {
        LAPLACE, ///< α = 2, μ = 0
        ASYMMETRIC_LAPLACE, ///< α = 2, μ ≠ 0
        CAUCHY, ///< α = 1, β = 0
        LEVY, ///< α = 0.5, |β| = 1
        UNITY_EXPONENT, ///< α = 1, β ≠ 0
        ONEHALF_EXPONENT, ///< α = 0.5, |β| ≠ 1
        GENERAL ///< the rest
    };

    DISTRIBUTION_TYPE distributionType = LAPLACE; ///< type of distribution (Laplace by default)

public:
    [[nodiscard]] inline double GetExponent() const noexcept { return alpha; }
    [[nodiscard]] inline double GetSkewness() const noexcept { return beta; }
    [[nodiscard]] inline double GetScale() const noexcept { return gamma; }
    [[nodiscard]] inline double GetLogScale() const noexcept { return logGamma; }

    [[nodiscard]] Support SupportType() const override;
    [[nodiscard]] RealType MinValue() const override;
    [[nodiscard]] RealType MaxValue() const override;

protected:
    double pdfLaplace(double x) const;
    double logpdfLaplace(double x) const;
    double cdfLaplace(double x) const;
    double cdfLaplaceCompl(double x) const;

private:
    double pdfByLevy(double x) const;
    double pdfByCauchy(double x) const;

    /**
     * @fn conditionalShift
     * @param w value of the exponential mixing variable
     * @return location of the conditional stable law of Y given W = w
     */
    double conditionalShift(double w) const;
    /**
     * @fn lowerLogMixingVariable
     * @param x0 argument minus shift
     * @return log of the smallest mixing variable w that matters for the density and the
     * distribution function at x0
     */
    double lowerLogMixingVariable(double x0) const;
    /**
     * @fn mixtureIntegral
     * @param g function of the mixing variable w
     * @param logWMin log of the smallest w that matters
     * @return E[g(W)] for W ~ Exp(1)
     */
    template <typename Function>
    double mixtureIntegral(const Function &g, double logWMin) const;

public:
    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;

private:
    template <std::uniform_random_bit_generator G>
    double variateForUnityExponent(double z, G &g) const;
    template <std::uniform_random_bit_generator G>
    double variateForGeneralExponent(double z, G &g) const;
    template <std::uniform_random_bit_generator G>
    double variateForOneHalfExponent(double z, G &g) const;
    template <std::uniform_random_bit_generator G>
    double variateByCauchy(double z, G &g) const;

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
     * @fn quantileLaplace
     * @param p input parameter in the interval (0, 1)
     * @return quantile for (Asymmetric) Laplace distribution
     */
    RealType quantileLaplace(double p) const;
    /**
     * @fn quantileLaplace1m
     * @param p input parameter in the interval (0, 1)
     * @return quantile of 1-p for (Asymmetric) Laplace distribution
     */
    RealType quantileLaplace1m(double p) const;

    /**
     * @fn searchTailQuantile
     * Quantile search that starts far tails from their power-law asymptote
     */
    RealType searchTailQuantile(double p, bool upperTail) const;
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;
};

/**
 * @brief The GeometricStableRand class <BR>
 * Geometric-Stable distribution
 *
 * Notation: X ~ GS(α, β, γ, μ)
 *
 * If X ~ Laplace(m, γ, κ), then X - m ~ GS(2, β, γ, (1/κ - κ) * γ) with arbitrary β
 */
template <typename RealType>
class GeometricStableRand final : public GeneralGeometricStableDistribution<RealType>
{
public:
    GeometricStableRand(double exponent = 2, double skewness = 0, double scale = randlib::consts::sqrt2,
                        double location = 0);
    virtual ~GeometricStableRand() {}

    [[nodiscard]] std::string Name() const override;

private:
    void ChangeAsymmetry();

public:
    void SetParameters(double exponent, double skewness);
    void SetLocation(double location);
    void SetScale(double scale);
    [[nodiscard]] inline double GetLocation() const noexcept { return this->mu; }
};

} // namespace randlib

#include <randlib/continuous/laplace.hpp>
#include <randlib/continuous/levy.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/cauchy.hpp>

namespace randlib
{

template <typename RealType>
GeneralGeometricStableDistribution<RealType>::GeneralGeometricStableDistribution(double exponent, double skewness,
                                                                                 double scale, double location,
                                                                                 double shift)
{
    SetParameters(exponent, skewness, scale, location, shift);
}

template <typename RealType>
void GeneralGeometricStableDistribution<RealType>::SetParameters(double exponent, double skewness, double scale,
                                                                 double location, double shift)
{
    if (!(exponent >= 0.1 && exponent <= 2.0))
        throw std::invalid_argument(
            "Geometric-Stable distribution: exponent should be inside the interval [0.1, 2], but it's equal to " +
            std::to_string(exponent));
    if (!(std::fabs(skewness) <= 1.0))
        throw std::invalid_argument(
            "Geometric-Stable distribution: skewness should be inside the interval [-1, 1], but it's equal to " +
            std::to_string(skewness));

    alpha = exponent;
    alphaInv = 1.0 / alpha;
    beta = skewness;
    Z.SetExponent(alpha);
    Z.SetSkewness(beta);

    SetScale(scale);
    SetLocation(location);
    SetShift(shift);

    if (alpha == 2.0)
        distributionType = (mu == 0.0) ? LAPLACE : ASYMMETRIC_LAPLACE;
    else if (alpha == 0.5)
        distributionType = (std::fabs(beta) == 1.0) ? LEVY : ONEHALF_EXPONENT;
    else if (alpha == 1.0)
        distributionType = (beta == 0.0) ? CAUCHY : UNITY_EXPONENT;
    else
        distributionType = GENERAL;
}

template <typename RealType>
void GeneralGeometricStableDistribution<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Geometric-Stable distribution: location should be finite");
    mu = location;
}

template <typename RealType>
void GeneralGeometricStableDistribution<RealType>::SetShift(double shift)
{
    if (!std::isfinite(shift))
        throw std::invalid_argument("Geometric-Stable distribution: shift should be finite");
    m = shift;
}

template <typename RealType>
void GeneralGeometricStableDistribution<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Asymmetric-Laplace distribution: scale should be positive, but it's equal to " +
                                    std::to_string(scale));
    gamma = scale;
    logGamma = std::log(gamma);
}

template <typename RealType>
void GeneralGeometricStableDistribution<RealType>::SetAsymmetry(double asymmetry)
{
    if (!(asymmetry > 0.0))
        throw std::invalid_argument(
            "Asymmetric-Laplace distribution: asymmetry parameter should be positive, but it's equal to " +
            std::to_string(asymmetry));
    kappa = asymmetry;
    kappaInv = 1.0 / kappa;
    kappaSq = kappa * kappa;
    log1pKappaSq = std::log1pl(kappaSq);
    logKappa = std::log(kappa);
}

template <typename RealType>
Support GeneralGeometricStableDistribution<RealType>::SupportType() const
{
    if (alpha < 1) {
        if (beta == 1 && mu >= 0)
            return Support::RightSemiFinite;
        if (beta == -1 && mu <= 0)
            return Support::LeftSemiFinite;
    }
    return Support::Infinite;
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::MinValue() const
{
    if (alpha < 1 && beta == 1 && mu >= 0)
        return m;
    return -INFINITY;
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::MaxValue() const
{
    if (alpha < 1 && beta == -1 && mu <= 0)
        return m;
    return INFINITY;
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::pdfLaplace(double x) const
{
    return std::exp(logpdfLaplace(x));
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::logpdfLaplace(double x) const
{
    double y = x / gamma;
    y *= (x < 0) ? kappaInv : -kappa;
    return y - logGamma - log1pKappaSq + logKappa;
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::cdfLaplace(double x) const
{
    double y = x / gamma;
    if (x < 0) {
        y *= kappaInv;
        y += 2 * logKappa - log1pKappaSq;
        return std::exp(y);
    }
    return -std::expm1l(-log1pKappaSq - kappa * y);
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::cdfLaplaceCompl(double x) const
{
    double y = x / gamma;
    if (x < 0) {
        y *= kappaInv;
        y += 2 * logKappa - log1pKappaSq;
        return -std::expm1l(y);
    }
    return std::exp(-log1pKappaSq - kappa * y);
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::pdfByLevy(double x) const
{
    if (mu == 0) {
        /// Invert parameter for case of -Levy
        x *= beta;
        if (x < 0)
            return 0.0;
        if (x == 0)
            return INFINITY;
        double halfSigmaInv = 0.5 / gamma;
        double z = std::sqrt(halfSigmaInv * x);
        double y = randlib::consts::inv_sqrtpi - math::xexpxsqerfc(z);
        return 0.5 * y / (gamma * z);
    }

    if (x == 0) {
        return std::fabs(1.0 / mu);
    }

    double a = x / gamma, c = mu / gamma;
    if (beta < 0) {
        a = -a;
        c = -c;
    }

    if (a < 0) {
        if (c >= 0.0)
            return 0.0;
        double h = std::sqrt(1 - 2 * c); // can be hashed
        double hm1 = h - 1;
        double y = 2 * a / (hm1 * hm1);
        y -= std::log(gamma * hm1 * h);
        y += randlib::consts::ln2;
        return std::exp(y);
    }

    double sqrt2a = std::sqrt(2 * a);

    if (c < 0.5) {
        double h = std::sqrt(1 - 2 * c); // can be hashed
        double term1 = math::sign(c) * sqrt2a / (1 - h);
        double term2 = sqrt2a / (1 + h);
        double y = math::xexpxsqerfc(term1);
        y -= math::xexpxsqerfc(term2);
        y /= h * gamma * sqrt2a;
        return y;
    }

    if (c == 0.5) {
        double y = math::xexpxsqerfc(sqrt2a);
        y *= (8.0 * a + 2.0) / sqrt2a;
        y -= 4 * sqrt2a / randlib::consts::sqrtpi;
        return y / gamma;
    }

    /// we do numerical integration in the case of c > 0.5
    // (as I wasn't able to handle complex integrals)
    double adivc = x / mu;
    double y = math::integrate(
        [a, c, adivc](double t) {
            if (t <= 0 || t >= adivc)
                return 0.0;
            double amct = a - c * t;
            double integrand = 0.5 * t / amct;
            ++integrand;
            integrand *= t;
            integrand = std::exp(-integrand);
            integrand *= t;
            integrand /= amct * std::sqrt(amct);
            return integrand;
        },
        0, adivc);
    return y / (gamma * randlib::consts::sqrt2pi);
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::pdfByCauchy(double x) const
{
    // TODO: find analytical solution, maybe using residue technique
    // otherwise, use peak, which is t = exp(x / sqrt(mu^2 + sigma^2))
    // to split the integral on two
    if (x == 0)
        return INFINITY;
    return -gamma * randlib::consts::inv_pi *
           math::integrate(
               [this, x](double t) {
                   if (t <= 0.0 || t >= 1.0)
                       return 0.0;
                   double logt = std::log(t);
                   double a = gamma * logt;
                   a *= a;
                   double b = x + mu * logt;
                   b *= b;
                   return logt / (a + b);
               },
               0, 1);
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::conditionalShift(double w) const
{
    /// Given W = w, Y = S(w) ~ S1(α, β, γ w^(1/α), μ w); for α = 1 the S1 scale family also
    /// shifts the location by (2/π) β γ w log(γ w)
    double shift = mu * w;
    if (distributionType == UNITY_EXPONENT)
        shift += randlib::consts::two_inv_pi * beta * gamma * w * std::log(gamma * w);
    return shift;
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::lowerLogMixingVariable(double x0) const
{
    /// Given W = w, (x0 - shift) / (γ w^(1/α)) leaves the bulk of the stable law below
    /// w = (|x0| / γ)^α; twenty e-folds below that, the integrands no longer change
    if (x0 == 0)
        return -40.0;
    return std::min(-20.0, alpha * std::log(std::fabs(x0) / gamma) - 20.0);
}

template <typename RealType>
template <typename Function>
double GeneralGeometricStableDistribution<RealType>::mixtureIntegral(const Function &g, double logWMin) const
{
    /// ∫ e^(-w) g(w) dw over W ~ Exp(1), integrated in u = log(w): the features of g have unit
    /// width in u even in the far tails, where they sit at w ~ |x|^α. Below w = exp(logWMin),
    /// g is taken constant; above w = e^7, e^(-w) underflows
    static constexpr double LOG_W_MAX = 7.0;
    double wMin = std::exp(logWMin);
    double head = -std::expm1(-wMin) * g(wMin);
    return head + math::integrate(
                      [&g](double u) {
                          double w = std::exp(u);
                          return w * std::exp(-w) * g(w);
                      },
                      logWMin, LOG_W_MAX, 1e-13);
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::f(const RealType &x) const
{
    double x0 = x - m;
    /// Cut zeros for semi-infinite distribution
    if (alpha < 1 && std::fabs(beta) == 1) {
        if (x0 < 0 && mu >= 0 && beta > 0)
            return 0.0;
        if (x0 > 0 && mu <= 0 && beta < 0)
            return 0.0;
    }

    /// Verify if any of analytical expressions is applicable
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE)
        return pdfLaplace(x0);
    if (distributionType == CAUCHY)
        return pdfByCauchy(x0);
    if (distributionType == LEVY)
        return pdfByLevy(x0);

    if (x0 == 0) {
        if (alpha == 1)
            return INFINITY;
        if (mu == 0) {
            double coef = std::tgammal(1.0 - alphaInv);
            if (!std::isfinite(coef))
                return INFINITY;
            return Z.f(0) * coef / gamma;
        }
    }

    if (distributionType != UNITY_EXPONENT) {
        return mixtureIntegral(
            [this, x0](double w) {
                double tau = 1.0 / (gamma * std::pow(w, alphaInv));
                return tau * Z.f(tau * x0 - tau * conditionalShift(w));
            },
            lowerLogMixingVariable(x0));
    }

    /// For α = 1, small w would need the stable law far in its tails, where it is slow (#1);
    /// integrate in z = exp(-w) instead
    return math::integrate(
        [this, x0](double z) {
            if (z <= 0 || z >= 1)
                return 0.0;
            double logz = std::log(z);
            double tau = 1.0 / (gamma * std::pow(-logz, alphaInv));
            double par = tau * (x0 - conditionalShift(-logz));
            return tau * Z.f(par);
        },
        0, 1, 1e-13);
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::logf(const RealType &x) const
{
    return (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE) ? logpdfLaplace(x - m)
                                                                                   : std::log(f(x));
}

template <typename RealType>
double GeneralGeometricStableDistribution<RealType>::F(const RealType &x) const
{
    double x0 = x - m;
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE)
        return cdfLaplace(x0);
    if (distributionType != UNITY_EXPONENT) {
        return mixtureIntegral(
            [this, x0](double w) { return Z.F((x0 - conditionalShift(w)) / (gamma * std::pow(w, alphaInv))); },
            lowerLogMixingVariable(x0));
    }

    /// For α = 1, integrate in t = z / (1 - z) instead (see f)
    return math::integrate(
        [this, x0](double z) {
            if (z <= 0)
                return (x0 == 0.0) ? Z.F(-mu / gamma) : ((x0 < 0.0) ? 0.0 : 1.0);
            if (z >= 1)
                return Z.F(-mu / gamma);
            double denominator = 1.0 - z;
            double t = z / denominator;
            double par = x0 - conditionalShift(t);
            par /= std::pow(t, alphaInv) * gamma;
            double y = std::exp(-t) * Z.F(par);
            return y / (denominator * denominator);
        },
        0, 1, 1e-13);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
double GeneralGeometricStableDistribution<RealType>::variateForUnityExponent(double z, G &g) const
{
    double W = detail::standardExponential<RealType>(g);
    double X = z + 2 * beta * std::log(gamma * W) / randlib::consts::pi;
    X *= gamma;
    X += mu;
    return X * W;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
double GeneralGeometricStableDistribution<RealType>::variateForGeneralExponent(double z, G &g) const
{
    double W = detail::standardExponential<RealType>(g);
    double X = std::pow(W, alphaInv) * gamma * z;
    return mu * W + X;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
double GeneralGeometricStableDistribution<RealType>::variateForOneHalfExponent(double z, G &g) const
{
    double W = detail::standardExponential<RealType>(g);
    double X = mu + gamma * W * z;
    return X * W;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
double GeneralGeometricStableDistribution<RealType>::variateByCauchy(double z, G &g) const
{
    double W = detail::standardExponential<RealType>(g);
    double X = mu + gamma * z;
    return X * W;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType GeneralGeometricStableDistribution<RealType>::operator()(G &g) const
{
    switch (distributionType) {
    case LAPLACE:
        return gamma * LaplaceRand<RealType>::StandardVariate(g);
    case ASYMMETRIC_LAPLACE:
        return gamma * AsymmetricLaplaceRand<RealType>::StandardVariate(kappa, g);
    case ONEHALF_EXPONENT:
    case LEVY:
        return variateForOneHalfExponent(Z(g), g);
    case CAUCHY:
        return variateByCauchy(Z(g), g);
    case UNITY_EXPONENT:
        return variateForUnityExponent(Z(g), g);
    case GENERAL:
        return variateForGeneralExponent(Z(g), g);
    default:
        throw std::runtime_error("Shifted Geometric Stable distribution: invalid distribution type");
    }
}

template <typename RealType>
void GeneralGeometricStableDistribution<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    switch (distributionType) {
    case LAPLACE: {
        for (RealType &var : outputData)
            var = gamma * LaplaceRand<RealType>::StandardVariate(g);
    } break;
    case ASYMMETRIC_LAPLACE: {
        for (RealType &var : outputData)
            var = gamma * AsymmetricLaplaceRand<RealType>::StandardVariate(kappa, g);
    } break;
    case ONEHALF_EXPONENT:
    case LEVY: {
        Z.Sample(g, outputData);
        for (RealType &var : outputData)
            var = variateForOneHalfExponent(var, g);
    } break;
    case CAUCHY: {
        Z.Sample(g, outputData);
        for (RealType &var : outputData)
            var = variateByCauchy(var, g);
    } break;
    case UNITY_EXPONENT: {
        Z.Sample(g, outputData);
        for (RealType &var : outputData)
            var = variateForUnityExponent(var, g);
    } break;
    case GENERAL: {
        Z.Sample(g, outputData);
        for (RealType &var : outputData)
            var = variateForGeneralExponent(var, g);
    } break;
    default:
        throw std::runtime_error("Shifted Geometric Stable distribution: invalid distribution type");
    }
}

template <typename RealType>
long double GeneralGeometricStableDistribution<RealType>::Mean() const
{
    if (alpha > 1)
        return m + mu;
    if (beta == 1)
        return INFINITY;
    return (beta == -1) ? -INFINITY : NAN;
}

template <typename RealType>
long double GeneralGeometricStableDistribution<RealType>::Variance() const
{
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE)
        return mu * mu + 2 * gamma * gamma;
    return INFINITY;
}

template <typename RealType>
std::complex<double> GeneralGeometricStableDistribution<RealType>::CFImpl(double t) const
{
    double x = 0;
    if (alpha != 2 && beta != 0) {
        x = (alpha == 1) ? randlib::consts::two_inv_pi * std::log(t) : -std::tan(randlib::consts::half_pi * alpha);
        x *= beta;
    }
    double re = std::pow(gamma * t, alpha);
    std::complex<double> psi = std::complex<double>(re, re * x - mu * t);
    return 1.0 / (1.0 + psi);
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::Median() const
{
    if (distributionType == LAPLACE || distributionType == ASYMMETRIC_LAPLACE)
        return quantileLaplace(0.5);
    return ContinuousDistribution<RealType>::Median();
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::Mode() const
{
    // TODO: calculate mode for more cases
    if (alpha == 1 || alpha == 2)
        return m;
    return ContinuousDistribution<RealType>::Mode();
}

template <typename RealType>
long double GeneralGeometricStableDistribution<RealType>::Skewness() const
{
    if (distributionType != LAPLACE && distributionType != ASYMMETRIC_LAPLACE)
        return NAN;
    long double k4 = kappaSq * kappaSq;
    long double k6 = k4 * kappaSq;
    long double z = (k4 + 1);
    z *= std::sqrt(z);
    long double y = 2 * (1 - k6);
    return y / z;
}

template <typename RealType>
long double GeneralGeometricStableDistribution<RealType>::ExcessKurtosis() const
{
    if (distributionType != LAPLACE && distributionType != ASYMMETRIC_LAPLACE)
        return NAN;
    long double denominator = kappaSq + kappaInv * kappaInv;
    denominator *= denominator;
    return 6.0 - 12.0 / denominator;
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::searchTailQuantile(double p, bool upperTail) const
{
    double tail = std::min(p, 1.0 - p);
    double guess = (tail < 0.05) ? stableTailQuantile(tail, upperTail != (p > 0.5), alpha, beta, gamma, mu) : NAN;
    return this->searchQuantile(p, upperTail, guess, 0.05);
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::quantileImpl(double p) const
{
    return searchTailQuantile(p, false);
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::quantileImpl1m(double p) const
{
    return searchTailQuantile(p, true);
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::quantileLaplace(double p) const
{
    if (p < kappaSq / (1 + kappaSq)) {
        RealType q = p / kappaSq + p;
        q = std::log(q);
        q *= kappa * gamma;
        return m + q;
    } else {
        RealType q = (kappaSq + 1) * (1 - p);
        q = std::log(q);
        q *= gamma / kappa;
        return m - q;
    }
}

template <typename RealType>
RealType GeneralGeometricStableDistribution<RealType>::quantileLaplace1m(double p) const
{
    if (p > 1.0 / (1 + kappaSq)) {
        RealType pm1 = p - 1.0;
        RealType q = -pm1 / kappaSq - pm1;
        q = std::log(q);
        q *= kappa * gamma;
        return m + q;
    } else {
        RealType q = (kappaSq + 1) * p;
        q = std::log(q);
        q *= gamma / kappa;
        return m - q;
    }
}
template <typename RealType>
GeometricStableRand<RealType>::GeometricStableRand(double exponent, double skewness, double scale, double location)
    : GeneralGeometricStableDistribution<RealType>(exponent, skewness, scale, location)
{
    ChangeAsymmetry();
}

template <typename RealType>
std::string GeometricStableRand<RealType>::Name() const
{
    return "Geometric Stable(" + this->toStringWithPrecision(this->GetExponent()) + ", " +
           this->toStringWithPrecision(this->GetSkewness()) + ", " + this->toStringWithPrecision(this->GetScale()) +
           ", " + this->toStringWithPrecision(this->GetLocation()) + ")";
}

template <typename RealType>
void GeometricStableRand<RealType>::ChangeAsymmetry()
{
    if (this->distributionType != this->LAPLACE && this->distributionType != this->ASYMMETRIC_LAPLACE)
        return;
    double gamma2 = 2 * this->gamma;
    double asymmetry = this->mu * this->mu + gamma2 * gamma2;
    asymmetry = std::sqrt(asymmetry) - this->mu;
    asymmetry /= gamma2;
    this->SetAsymmetry(asymmetry);
}

template <typename RealType>
void GeometricStableRand<RealType>::SetParameters(double exponent, double skewness)
{
    GeneralGeometricStableDistribution<RealType>::SetParameters(exponent, skewness);
    ChangeAsymmetry();
}

template <typename RealType>
void GeometricStableRand<RealType>::SetLocation(double location)
{
    GeneralGeometricStableDistribution<RealType>::SetLocation(location);
    ChangeAsymmetry();
}

template <typename RealType>
void GeometricStableRand<RealType>::SetScale(double scale)
{
    GeneralGeometricStableDistribution<RealType>::SetScale(scale);
    ChangeAsymmetry();
}

} // namespace randlib

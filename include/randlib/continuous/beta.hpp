#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/core/exponential_family.hpp>

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN

namespace randlib
{

/**
 * @brief The BetaDistribution class <BR>
 * Abstract class for Beta distribution
 *
 * f(x | α, β) = y^{α-1} (1-y)^{β-1} / B(α, β), <BR>
 * where y = (x-a)/(b-a) and B(α, β) denotes Beta function
 *
 * Notation: X ~ B(α, β, a, b) or X ~ B(α, β) for a=0 and b=1
 *
 * Related distributions (a=0, b=1): <BR>
 * 1 − X ~ B(β, α) <BR>
 * X / (1 - X) ~ B'(α, β) <BR>
 * X = Y / (Y + Z), where Y ~ Γ(α) and Z ~ Γ(β) <BR>
 * βX / α(1 - X) ~ F(2α, 2β)
 */
template <typename RealType>
class BetaDistribution : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

protected:
    double alpha = 1; ///< first shape α
    double beta = 1; ///< second shape β
    double a = 0; ///< min bound
    double b = 1; ///< max bound
    double bma = 1; ///< b-a
    double bmaInv = 1; ///< 1/(b-a)
    double logbma = 0; ///< log(b-a)

    GammaRand<RealType> GammaRV1{}, GammaRV2{};

private:
    static constexpr double edgeForGenerators = 8.0;
    double logBetaFun = 0; ///< log(B(α, β)
    double betaFun = 1; ///< B(α, β)

    static double logStandardized(double distance, double complement);

    /// coefficients for generators
    struct genCoef_t
    {
        double s, t, u;
    } genCoef = {0, 0, 0};

protected:
    BetaDistribution(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1);

public:
    virtual ~BetaDistribution() {}

    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return a; }
    [[nodiscard]] RealType MaxValue() const override { return b; }

private:
    enum GENERATOR_ID
    {
        UNIFORM, ///< standard uniform variate
        ARCSINE, ///< arcsine method
        CHENG, ///< Cheng's method
        REJECTION_UNIFORM, ///< rejection method from uniform distribution for specific value of shapes α = β = 1.5
        REJECTION_UNIFORM_EXTENDED, ///< rejection method from uniform distribution accelerated by using exponential distribution
        REJECTION_NORMAL, ///< rejection method normal distribution
        JOHNK, ///< Johnk's method
        ATKINSON_WHITTAKER, ///< Atkinson-Whittaker's method
        GAMMA_RATIO ///< ratio of two gamma variables
    };

    /**
     * @return id of used variate generator according to the shapes
     */
    GENERATOR_ID getIdOfUsedGenerator() const
    {
        if (alpha < 1 && beta < 1 && alpha + beta > 1)
            return ATKINSON_WHITTAKER;

        if (math::detail::areClose(alpha, beta)) {
            if (math::detail::areClose(alpha, 1.0))
                return UNIFORM;
            else if (math::detail::areClose(alpha, 0.5))
                return ARCSINE;
            else if (math::detail::areClose(alpha, 1.5))
                return REJECTION_UNIFORM;
            else if (alpha > 1)
                return (alpha < 2) ? REJECTION_UNIFORM_EXTENDED : REJECTION_NORMAL;
        }
        if (std::min(alpha, beta) > 0.5 && std::max(alpha, beta) > 1)
            return CHENG;
        return (alpha + beta < 2) ? JOHNK : GAMMA_RATIO;
    }

    /**
     */
    void setCoefficientsForGenerator();

protected:
    /**
     * @param shape1 α
     * @param shape2 β
     */
    void SetShapes(double shape1, double shape2);
    /**
     * @param minValue a
     * @param maxValue b
     */
    void SetSupport(double minValue, double maxValue);

public:
    /**
     * @return first shape α
     */
    [[nodiscard]] inline double GetAlpha() const noexcept { return alpha; }
    /**
     * @return second shape β
     */
    [[nodiscard]] inline double GetBeta() const noexcept { return beta; }
    /**
     * @return B(α, β)
     */
    [[nodiscard]] inline double GetBetaFunction() const noexcept { return betaFun; }
    /**
     * @return log(B(α, β))
     */
    [[nodiscard]] inline double GetLogBetaFunction() const noexcept { return logBetaFun; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

private:
    /**
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for α = β = 1.5
     */
    template <std::uniform_random_bit_generator G>
    RealType variateRejectionUniform(G &g) const;

    /**
     * Symmetric beta generator via rejection from the uniform density
     * @return beta variate for 1 < α = β < 2 and α != 1.5
     */
    template <std::uniform_random_bit_generator G>
    RealType variateRejectionUniformExtended(G &g) const;

    /**
     * Arcsine beta generator
     * @return beta variate for α = β = 0.5
     */
    template <std::uniform_random_bit_generator G>
    RealType variateArcsine(G &g) const;

    /**
     * Symmetric beta generator via rejection from the normal density
     * @return beta variate for equal shape parameters > 2
     */
    template <std::uniform_random_bit_generator G>
    RealType variateRejectionNormal(G &g) const;

    /**
     * Johnk's beta generator
     * @return beta variate for small shape parameters < 1
     */
    template <std::uniform_random_bit_generator G>
    RealType variateJohnk(G &g) const;

    /**
     * Cheng's beta generator
     * @return beta variate for max(α, β) > 1 and min(α, β) > 0.5
     */
    template <std::uniform_random_bit_generator G>
    RealType variateCheng(G &g) const;

    /**
     * Atkinson-Whittaker beta generator
     * @return beta variate for max(α, β) < 1 and α + β > 1
     */
    template <std::uniform_random_bit_generator G>
    RealType variateAtkinsonWhittaker(G &g) const;

    /**
     * Gamma ratio beta generator
     * @return beta variate for the rest variations of shapes
     */
    template <std::uniform_random_bit_generator G>
    RealType variateGammaRatio(G &g) const;

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    /**
     * @return E[ln(X)]
     */
    long double GeometricMean() const;
    [[nodiscard]] long double Variance() const override;
    /**
     * @return Var(ln(X))
     */
    long double GeometricVariance() const;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
    /**
     * @brief MeanAbsoluteDeviation
     * @return E[|X - E[X]|]
     */
    long double MeanAbsoluteDeviation() const;

protected:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

    static constexpr char ALPHA_ZERO[] =
        "Possibly one or more elements of the sample coincide with the lower boundary a.";
    static constexpr char BETA_ZERO[] =
        "Possibly one or more elements of the sample coincide with the upper boundary b.";
};

/**
 * @brief The BetaRand class <BR>
 * Beta distribution
 */
template <typename RealType>
class BetaRand final : public BetaDistribution<RealType>, public ExponentialFamily<RealType, DoublePair>
{
public:
    BetaRand(double shape1 = 1, double shape2 = 1, double minValue = 0, double maxValue = 1)
        : BetaDistribution<RealType>(shape1, shape2, minValue, maxValue)
    {
    }
    [[nodiscard]] std::string Name() const override;

    using BetaDistribution<RealType>::SetShapes;
    using BetaDistribution<RealType>::SetSupport;

    DoublePair SufficientStatistic(RealType x) const override;
    DoublePair SourceParameters() const override;
    DoublePair SourceToNatural(DoublePair sourceParameters) const override;
    double LogNormalizer(DoublePair theta) const override;
    DoublePair LogNormalizerGradient(DoublePair theta) const override;
    double CarrierMeasure(RealType) const override;

    /**
     * @brief GetSampleLogMeanNorm
     * @param sample
     * @return mean average of log(x) for x from normalized sample
     */
    [[nodiscard]] long double GetSampleLogMeanNorm(std::span<const RealType> sample) const;
    /**
     * @brief GetSampleLog1pMeanNorm
     * @param sample
     * @return mean average of log(1+x) for x from normalized sample
     */
    [[nodiscard]] long double GetSampleLog1pMeanNorm(std::span<const RealType> sample) const;
    /**
     * @brief GetSampleLog1mMeanNorm
     * @param sample
     * @return mean average of log(1-x) for x from normalized sample
     */
    [[nodiscard]] long double GetSampleLog1mMeanNorm(std::span<const RealType> sample) const;

    /**
     * set α, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG normalized sample average of ln(X)
     * @param meanNorm normalized sample average
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitAlpha(long double lnG, long double meanNorm);

    /**
     * set α, estimated via maximum likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitAlpha(std::span<const RealType> sample);

    /**
     * set β, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG1m normalized sample average of ln(1-X)
     * @param meanNorm normalized sample average
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitBeta(long double lnG1m, long double meanNorm);

    /**
     * set β, estimated via maximum likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitBeta(std::span<const RealType> sample);

    /**
     * set α and β, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG sample average of ln(X)
     * @param lnG1m sample average of ln(1-X)
     * @param mean sample average
     * @param variance sample variance
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitShapes(long double lnG, long double lnG1m, long double mean,
                                                 long double variance);

    /**
     * set α and β, estimated via maximum likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitShapes(std::span<const RealType> sample);
};

/**
 * @brief The ArcsineRand class <BR>
 * Arcsine distribution
 *
 * Notation: X ~ Arcsine(α, a, b)
 * 
 * Related distributions: <BR>
 * X ~ B(1 - α, α, a, b)
 */
template <typename RealType>
class ArcsineRand final : public BetaDistribution<RealType>
{
public:
    ArcsineRand(double shape = 0.5, double minValue = 0, double maxValue = 1)
        : BetaDistribution<RealType>(1.0 - shape, shape, minValue, maxValue)
    {
    }
    [[nodiscard]] std::string Name() const override;

    using BetaDistribution<RealType>::SetSupport;

    void SetShape(double shape);
    [[nodiscard]] inline double GetShape() const noexcept { return this->beta; }

    /**
     * set α and β, estimated via maximum likelihood,
     * using sufficient statistics instead of the whole sample
     * @param lnG average of all ln(X)
     * @param lnG1m average of all ln(1-X)
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitShape(long double lnG, long double lnG1m);

    /**
     * set α, estimated via maximum likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitShape(std::span<const RealType> sample);
};

/**
 * @brief The BaldingNicholsRand class <BR>
 * Balding-Nichols distribution
 *
 * Notation: X ~ Balding-Nichols(F, p)
 *
 * Related distributions: <BR>
 * X ~ B(p * F', (1 - p) * F', 0, 1) for F' = (1 - F) / F
 */
template <typename RealType>
class BaldingNicholsRand final : public BetaDistribution<RealType>
{
    double fixationIndex = 0.5, p = 0.5;

public:
    BaldingNicholsRand(double fixatingIndex, double frequency);
    [[nodiscard]] std::string Name() const override;

    void SetFixatingIndexAndFrequency(double fixatingIndex, double frequency);
    [[nodiscard]] inline double GetFixatingIndex() const noexcept { return fixationIndex; }
    [[nodiscard]] inline double GetFrequency() const noexcept { return p; }
};

} // namespace randlib

#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/discrete/bernoulli.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/exponential.hpp>

namespace randlib
{

template <typename RealType>
BetaDistribution<RealType>::BetaDistribution(double shape1, double shape2, double minValue, double maxValue)
{
    SetShapes(shape1, shape2);
    SetSupport(minValue, maxValue);
}

template <typename RealType>
void BetaDistribution<RealType>::setCoefficientsForGenerator()
{
    GENERATOR_ID id = getIdOfUsedGenerator();
    if (id == REJECTION_NORMAL) {
        double alpham1 = alpha - 1;
        genCoef.s = alpham1 * std::log1pl(0.5 / alpham1) - 0.5;
        genCoef.t = 1.0 / std::sqrt(8 * alpha - 4);
    } else if (id == CHENG) {
        genCoef.s = alpha + beta;
        genCoef.t = std::min(alpha, beta);
        if (genCoef.t > 1)
            genCoef.t = std::sqrt((2 * alpha * beta - genCoef.s) / (genCoef.s - 2));
        genCoef.u = alpha + genCoef.t;
    } else if (id == ATKINSON_WHITTAKER) {
        genCoef.t = std::sqrt(alpha * (1 - alpha));
        genCoef.t /= (genCoef.t + std::sqrt(beta * (1 - beta)));
        genCoef.s = beta * genCoef.t;
        genCoef.s /= (genCoef.s + alpha * (1 - genCoef.t));
    }
}

template <typename RealType>
void BetaDistribution<RealType>::SetShapes(double shape1, double shape2)
{
    if (!(shape1 > 0 && shape2 > 0))
        throw std::invalid_argument("Beta distribution: shapes should be positive");
    GammaRV1.SetParameters(shape1, 1);
    GammaRV2.SetParameters(shape2, 1);
    alpha = GammaRV1.GetShape();
    beta = GammaRV2.GetShape();
    logBetaFun = math::logBeta(alpha, beta);
    betaFun = std::exp(logBetaFun);
    setCoefficientsForGenerator();
}

template <typename RealType>
void BetaDistribution<RealType>::SetSupport(double minValue, double maxValue)
{
    if (!(minValue < maxValue))
        throw std::invalid_argument("Beta distribution: minimum value should be smaller than maximum value");

    a = minValue;
    b = maxValue;
    bma = b - a;
    bmaInv = 1.0 / bma;
    logbma = std::log(bma);
}

template <typename RealType>
double BetaDistribution<RealType>::f(const RealType &x) const
{
    if (x < a || x > b)
        return 0.0;
    if (x == a) {
        if (alpha == 1)
            return beta / bma;
        return (alpha > 1) ? 0 : INFINITY;
    }
    if (x == b) {
        if (beta == 1)
            return alpha / bma;
        return (beta > 1) ? 0 : INFINITY;
    }
    return std::exp(logf(x));
}

template <typename RealType>
double BetaDistribution<RealType>::logf(const RealType &x) const
{
    if (x < a || x > b)
        return -INFINITY;
    if (x == a) {
        if (alpha == 1)
            return std::log(beta / bma);
        return (alpha > 1) ? -INFINITY : INFINITY;
    }
    if (x == b) {
        if (beta == 1)
            return std::log(alpha / bma);
        return (beta > 1) ? -INFINITY : INFINITY;
    }
    const double xSt = (x - a) * bmaInv, ySt = (b - x) * bmaInv;
    double y = (alpha - 1) * logStandardized(xSt, ySt);
    y += (beta - 1) * logStandardized(ySt, xSt);
    return y - logBetaFun - logbma;
}

template <typename RealType>
double BetaDistribution<RealType>::logStandardized(double distance, double complement)
{
    /// log of a standardized distance, through its complement when that one is the accurate one
    return (complement < 0.5) ? std::log1p(-complement) : std::log(distance);
}

template <typename RealType>
double BetaDistribution<RealType>::F(const RealType &x) const
{
    if (x <= a)
        return 0.0;
    if (x >= b)
        return 1.0;
    /// Standardize: distances to both ends are exact, unlike 1 minus the standardized value
    double xSt = std::min((x - a) * bmaInv, 1.0), ySt = std::min((b - x) * bmaInv, 1.0);
    /// Workaround known case
    if (alpha == beta && beta == 0.5)
        return (xSt <= ySt) ? randlib::consts::two_inv_pi * std::asin(std::sqrt(xSt))
                            : 1.0 - randlib::consts::two_inv_pi * std::asin(std::sqrt(ySt));
    return math::ibeta(xSt, alpha, beta, logBetaFun, logStandardized(xSt, ySt), logStandardized(ySt, xSt));
}

template <typename RealType>
double BetaDistribution<RealType>::S(const RealType &x) const
{
    if (x <= a)
        return 1.0;
    if (x >= b)
        return 0.0;
    /// Standardize: distances to both ends are exact, unlike 1 minus the standardized value
    double xSt = std::min((x - a) * bmaInv, 1.0), ySt = std::min((b - x) * bmaInv, 1.0);
    /// Workaround known case
    if (alpha == beta && beta == 0.5)
        return (ySt <= xSt) ? randlib::consts::two_inv_pi * std::asin(std::sqrt(ySt))
                            : 1.0 - randlib::consts::two_inv_pi * std::asin(std::sqrt(xSt));
    return math::ibeta(ySt, beta, alpha, logBetaFun, logStandardized(ySt, xSt), logStandardized(xSt, ySt));
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateArcsine(G &g) const
{
    RealType U = 2 * detail::uniform01<RealType>(g) - 1;
    RealType X = std::sin(randlib::consts::pi * U);
    return X * X;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateRejectionUniform(G &g) const
{
    size_t iter = 0;
    do {
        RealType U = detail::uniform01<RealType>(g);
        RealType V = detail::uniform01<RealType>(g);
        if (0.25 * V * V <= U - U * U)
            return U;
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateRejectionUniformExtended(G &g) const
{
    size_t iter = 0;
    static constexpr double M_LN4 = randlib::consts::ln2 + randlib::consts::ln2;
    do {
        RealType U = detail::uniform01<RealType>(g);
        RealType W = detail::standardExponential<RealType>(g);
        RealType edge = M_LN4 + std::log(U - U * U);
        if (W >= (1.0 - alpha) * edge)
            return U;
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateRejectionNormal(G &g) const
{
    size_t iter = 0;
    RealType N = 0, Z = 0;
    RealType alpham1 = alpha - 1;
    RealType alpha2m1 = alpha + alpham1;
    do {
        do {
            N = detail::standardNormal<RealType>(g);
            Z = N * N;
        } while (Z >= alpha2m1);

        RealType W = detail::standardExponential<RealType>(g) + genCoef.s;
        RealType aux = 0.5 - alpham1 / (alpha2m1 - Z);
        aux *= Z;
        if (W + aux >= 0)
            return 0.5 + N * genCoef.t;
        aux = std::log1pl(-Z / alpha2m1);
        aux *= alpham1;
        aux += W + 0.5 * Z;
        if (aux >= 0)
            return 0.5 + N * genCoef.t;
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateJohnk(G &g) const
{
    RealType X = 0, Z = 0;
    RealType W = 0, V = 0;
    do {
        W = detail::standardExponential<RealType>(g) / alpha;
        V = detail::standardExponential<RealType>(g) / beta;
        X = std::exp(-W);
        Z = X + std::exp(-V);
    } while (Z > 1);
    return (Z > 0) ? (X / Z) : (W < V);
    ;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateCheng(G &g) const
{
    RealType R, T, Y;
    do {
        RealType U = detail::uniform01<RealType>(g);
        RealType V = detail::uniform01<RealType>(g);
        RealType X = std::log(U / (1 - U)) / genCoef.t;
        Y = alpha * std::exp(X);
        R = 1.0 / (beta + Y);
        T = 4 * U * U * V;
        T = std::log(T);
        T -= genCoef.u * X;
        T -= genCoef.s * std::log(genCoef.s * R);
    } while (T > 0);
    return Y * R;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateAtkinsonWhittaker(G &g) const
{
    size_t iter = 0;
    do {
        RealType U = detail::uniform01<RealType>(g);
        RealType W = detail::standardExponential<RealType>(g);
        if (U <= genCoef.s) {
            RealType X = genCoef.t * std::pow(U / genCoef.s, 1.0 / alpha);
            if (W >= (1.0 - beta) * std::log((1.0 - X) / (1.0 - genCoef.t)))
                return X;
        } else {
            RealType X = 1.0 - (1.0 - genCoef.t) * std::pow((1.0 - U) / (1.0 - genCoef.s), 1.0 / beta);
            if (W >= (1.0 - alpha) * std::log(X / genCoef.t))
                return X;
        }
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::variateGammaRatio(G &g) const
{
    RealType Y = GammaRV1(g);
    RealType Z = GammaRV2(g);
    return Y / (Y + Z);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaDistribution<RealType>::operator()(G &g) const
{
    double var = 0;
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM:
        var = detail::uniform01<RealType>(g);
        break;
    case ARCSINE:
        var = variateArcsine(g);
        break;
    case CHENG:
        var = variateCheng(g);
        break;
    case REJECTION_UNIFORM:
        var = variateRejectionUniform(g);
        break;
    case REJECTION_UNIFORM_EXTENDED:
        var = variateRejectionUniformExtended(g);
        break;
    case REJECTION_NORMAL:
        var = variateRejectionNormal(g);
        break;
    case JOHNK:
        var = variateJohnk(g);
        break;
    case ATKINSON_WHITTAKER:
        var = variateAtkinsonWhittaker(g);
        break;
    case GAMMA_RATIO:
    default:
        var = variateGammaRatio(g);
        break;
    }

    return a + bma * var;
}

template <typename RealType>
void BetaDistribution<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM: {
        for (RealType &var : outputData)
            var = detail::uniform01<RealType>(g);
    } break;
    case ARCSINE: {
        for (RealType &var : outputData)
            var = variateArcsine(g);
    } break;
    case CHENG: {
        for (RealType &var : outputData)
            var = variateCheng(g);
    } break;
    case REJECTION_UNIFORM: {
        for (RealType &var : outputData)
            var = variateRejectionUniform(g);
    } break;
    case REJECTION_UNIFORM_EXTENDED: {
        for (RealType &var : outputData)
            var = variateRejectionUniformExtended(g);
    } break;
    case REJECTION_NORMAL: {
        for (RealType &var : outputData)
            var = variateRejectionNormal(g);
    } break;
    case JOHNK: {
        for (RealType &var : outputData)
            var = variateJohnk(g);
    } break;
    case ATKINSON_WHITTAKER: {
        for (RealType &var : outputData)
            var = variateAtkinsonWhittaker(g);
    } break;
    case GAMMA_RATIO:
    default: {
        GammaRV1.Sample(g, outputData);
        for (RealType &var : outputData)
            var /= (var + GammaRV2(g));
    } break;
    }

    /// Shift and scale
    for (RealType &var : outputData)
        var = a + bma * var;
}

template <typename RealType>
long double BetaDistribution<RealType>::Mean() const
{
    double mean = alpha / (alpha + beta);
    return a + bma * mean;
}

template <typename RealType>
long double BetaDistribution<RealType>::GeometricMean() const
{
    return math::digamma(alpha) - math::digamma(alpha + beta);
}

template <typename RealType>
long double BetaDistribution<RealType>::Variance() const
{
    double var = alpha + beta;
    var *= var * (var + 1);
    var = alpha * beta / var;
    return bma * bma * var;
}

template <typename RealType>
long double BetaDistribution<RealType>::GeometricVariance() const
{
    return math::trigamma(alpha) - math::trigamma(alpha + beta);
}

template <typename RealType>
RealType BetaDistribution<RealType>::Median() const
{
    if (alpha == beta)
        return a + bma * 0.5;
    if (alpha == 1.0)
        return a - bma * std::expm1l(-randlib::consts::ln2 / beta);
    if (beta == 1.0)
        return a + bma * std::exp(-randlib::consts::ln2 / alpha);
    if (alpha >= 1.0 && beta >= 1.0) {
        double initValue = 3 * alpha - 1.0;
        initValue /= 3 * (alpha + beta) - 2.0;
        initValue *= bma;
        initValue += a;
        return ContinuousDistribution<RealType>::quantileImpl(0.5, initValue);
    }
    return ContinuousDistribution<RealType>::quantileImpl(0.5);
}

template <typename RealType>
RealType BetaDistribution<RealType>::Mode() const
{
    double mode;
    if (alpha > 1)
        mode = (beta > 1) ? (alpha - 1) / (alpha + beta - 2) : 1.0;
    else
        mode = (beta > 1) ? 0.0 : (alpha > beta);
    return a + bma * mode;
}

template <typename RealType>
long double BetaDistribution<RealType>::Skewness() const
{
    long double skewness = (alpha + beta + 1) / (alpha * beta);
    skewness = std::sqrt(skewness);
    skewness *= beta - alpha;
    skewness /= alpha + beta + 2;
    return 2 * skewness;
}

template <typename RealType>
long double BetaDistribution<RealType>::ExcessKurtosis() const
{
    long double sum = alpha + beta;
    long double kurtosis = alpha - beta;
    kurtosis *= kurtosis;
    kurtosis *= (sum + 1);
    kurtosis /= (alpha * beta * (sum + 2));
    --kurtosis;
    kurtosis /= (sum + 3);
    return 6 * kurtosis;
}

template <typename RealType>
long double BetaDistribution<RealType>::MeanAbsoluteDeviation() const
{
    double y = randlib::consts::ln2;
    y += alpha * std::log(alpha);
    y += beta * std::log(beta);
    y -= (alpha + beta + 1) * std::log(alpha + beta);
    y -= logBetaFun;
    y += logbma;
    return std::exp(y);
}

template <typename RealType>
RealType BetaDistribution<RealType>::quantileImpl(double p) const
{
    if (alpha == beta) {
        if (alpha == 0.5) {
            double x = std::sin(0.5 * randlib::consts::pi * p);
            return a + bma * x * x;
        }
        if (alpha == 1.0)
            return a + bma * p;
    }
    if (alpha == 1.0)
        return a - bma * std::expm1l(std::log1pl(-p) / beta);
    if (beta == 1.0)
        return a + bma * std::pow(p, 1.0 / alpha);
    return ContinuousDistribution<RealType>::quantileImpl(p);
}

template <typename RealType>
RealType BetaDistribution<RealType>::quantileImpl1m(double p) const
{
    if (alpha == beta) {
        if (alpha == 0.5) {
            double x = std::cos(0.5 * randlib::consts::pi * p);
            return a + bma * x * x;
        }
        if (alpha == 1.0)
            return b - bma * p;
    }
    if (alpha == 1.0)
        return a - bma * std::expm1l(std::log(p) / beta);
    if (beta == 1.0)
        return a + bma * std::exp(std::log1pl(-p) / alpha);
    return ContinuousDistribution<RealType>::quantileImpl1m(p);
}

template <typename RealType>
std::complex<double> BetaDistribution<RealType>::CFImpl(double t) const
{
    /// if we don't have singularity points, we can use direct integration
    if (alpha >= 1 && beta >= 1)
        return UnivariateDistribution<RealType>::CFImpl(t);

    double z = bma * t;
    double sinZ = std::sin(z);
    double cosZm1 = std::cos(z) - 1.0;

    double re = math::integrate(
        [this, z, cosZm1](double x) {
            if (x >= 1)
                return 0.0;
            if (x <= 0)
                return -cosZm1;
            double y = std::cos(z * x) - 1;
            y *= std::pow(x, alpha - 1);
            y -= cosZm1;
            return std::pow(1.0 - x, beta - 1) * y;
        },
        0, 1);
    re += betaFun;
    re += cosZm1 / beta;

    double im = math::integrate(
        [this, z, sinZ](double x) {
            if (x >= 1)
                return 0.0;
            if (x <= 0)
                return -sinZ;
            double y = std::sin(z * x);
            y *= std::pow(x, alpha - 1);
            y -= sinZ;
            return std::pow(1.0 - x, beta - 1) * y;
        },
        0, 1);
    im += sinZ / beta;

    std::complex<double> y(re, im);
    double cosTA = std::cos(t * a), sinTA = std::sin(t * a);
    return y * std::complex<double>(cosTA, sinTA) / betaFun;
}

template <typename RealType>
constexpr char BetaDistribution<RealType>::ALPHA_ZERO[];
template <typename RealType>
constexpr char BetaDistribution<RealType>::BETA_ZERO[];
// BETARAND

template <typename RealType>
std::string BetaRand<RealType>::Name() const
{
    return "Beta(" + this->toStringWithPrecision(this->GetAlpha()) + ", " +
           this->toStringWithPrecision(this->GetBeta()) + ", " + this->toStringWithPrecision(this->MinValue()) + ", " +
           this->toStringWithPrecision(this->MaxValue()) + ")";
}

template <typename RealType>
DoublePair BetaRand<RealType>::SufficientStatistic(RealType x) const
{
    double y = (x - this->a) * this->bmaInv;
    return {std::log(y), std::log1pl(-y)};
}

template <typename RealType>
DoublePair BetaRand<RealType>::SourceParameters() const
{
    return {this->alpha, this->beta};
}

template <typename RealType>
DoublePair BetaRand<RealType>::SourceToNatural(DoublePair sourceParameters) const
{
    return {sourceParameters.first - 1, sourceParameters.second - 1};
}

template <typename RealType>
double BetaRand<RealType>::LogNormalizer(DoublePair theta) const
{
    return this->logbma + math::logBeta(theta.first + 1, theta.second + 1);
}

template <typename RealType>
DoublePair BetaRand<RealType>::LogNormalizerGradient(DoublePair theta) const
{
    double psi1 = math::digamma(theta.first + 1);
    double psi2 = math::digamma(theta.second + 1);
    double psisum = math::digamma(theta.first + theta.second + 2);
    return {psi1 - psisum, psi2 - psisum};
}

template <typename RealType>
double BetaRand<RealType>::CarrierMeasure(RealType) const
{
    return 0;
}

template <typename RealType>
long double BetaRand<RealType>::GetSampleLogMeanNorm(std::span<const RealType> sample) const
{
    long double lnG = 0;
    for (RealType var : sample) {
        RealType x = (var - this->a) * this->bmaInv;
        lnG += std::log(x);
    }
    return lnG / sample.size();
}

template <typename RealType>
long double BetaRand<RealType>::GetSampleLog1pMeanNorm(std::span<const RealType> sample) const
{
    long double lnG1p = 0;
    for (RealType var : sample) {
        RealType x = (var - this->a) * this->bmaInv;
        lnG1p += std::log1pl(x);
    }
    return lnG1p / sample.size();
}

template <typename RealType>
long double BetaRand<RealType>::GetSampleLog1mMeanNorm(std::span<const RealType> sample) const
{
    long double lnG1m = 0;
    for (RealType var : sample) {
        RealType x = (var - this->a) * this->bmaInv;
        lnG1m += std::log1pl(-x);
    }
    return lnG1m / sample.size();
}

template <typename RealType>
void BetaRand<RealType>::FitAlpha(long double lnG, long double meanNorm)
{
    if (meanNorm <= 0 || meanNorm >= 1)
        throw this->fitError(FitErrorCode::NotApplicable,
                             "Normalized mean of the sample should be in interval of (0, 1)");
    if (this->beta == 1.0) {
        /// for β = 1 we have explicit expression for estimator
        SetShapes(-1.0 / lnG, this->beta);
    } else {
        /// get initial value for shape by method of moments
        RealType shape = meanNorm;
        shape /= 1.0 - meanNorm;
        shape *= this->beta;
        /// run root-finding procedure
        if (!math::findRootNewtonFirstOrder<RealType>(
                [this, lnG](double x) {
                    double first = math::digamma(x) - math::digamma(x + this->beta) - lnG;
                    double second = math::trigamma(x);
                    return DoublePair(first, second);
                },
                shape))
            throw this->fitError(FitErrorCode::Undefined, "Error in root-finding procedure.");
        SetShapes(shape, this->beta);
    }
}

template <typename RealType>
void BetaRand<RealType>::FitAlpha(std::span<const RealType> sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b));

    long double lnG = this->GetSampleLogMeanNorm(sample);
    if (!std::isfinite(lnG))
        throw this->fitError(FitErrorCode::WrongReturn, this->ALPHA_ZERO);
    long double mean = this->GetSampleMean(sample);
    mean -= this->a;
    mean *= this->bmaInv;
    FitAlpha(lnG, mean);
}

template <typename RealType>
void BetaRand<RealType>::FitBeta(long double lnG1m, long double meanNorm)
{
    if (meanNorm <= 0 || meanNorm >= 1)
        throw this->fitError(FitErrorCode::NotApplicable,
                             "Normalized mean of the sample should be in interval of (0, 1)");
    if (this->alpha == 1.0) {
        /// for α = 1 we have explicit expression for estimator
        SetShapes(this->alpha, -1.0 / lnG1m);
    } else {
        /// get initial value for shape by method of moments
        RealType shape = this->alpha / meanNorm - this->alpha;
        /// run root-finding procedure
        if (!math::findRootNewtonFirstOrder<RealType>(
                [this, lnG1m](double x) {
                    double first = math::digamma(x) - math::digamma(x + this->alpha) - lnG1m;
                    double second = math::trigamma(x);
                    return DoublePair(first, second);
                },
                shape))
            throw this->fitError(FitErrorCode::Undefined, "Error in root-finding procedure.");
        this->SetShapes(this->alpha, shape);
    }
}

template <typename RealType>
void BetaRand<RealType>::FitBeta(std::span<const RealType> sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b));

    long double lnG1m = this->GetSampleLog1mMeanNorm(sample);
    if (!std::isfinite(lnG1m))
        throw this->fitError(FitErrorCode::WrongReturn, this->BETA_ZERO);
    long double mean = this->GetSampleMean(sample);
    mean -= this->a;
    mean *= this->bmaInv;
    FitBeta(lnG1m, mean);
}

template <typename RealType>
void BetaRand<RealType>::FitShapes(long double lnG, long double lnG1m, long double mean, long double variance)
{
    /// get initial values for shapes by method of moments
    double scaledMean = (mean - this->a) * this->bmaInv;
    double scaledVar = variance * this->bmaInv * this->bmaInv;
    double temp = scaledMean * (1.0 - scaledMean) / scaledVar - 1.0;
    double shape1 = 0.001, shape2 = shape1;
    if (temp > 0) {
        shape1 = scaledMean * temp;
        shape2 = (1.0 - scaledMean) * temp;
    }
    DoublePair shapes = std::make_pair(shape1, shape2);

    /// run root-finding procedure
    if (!math::findRootNewtonFirstOrder2d(
            [lnG, lnG1m](DoublePair x) {
                double digammaAlphapBeta = math::digamma(x.first + x.second);
                double digammaAlpha = math::digamma(x.first);
                double digammaBeta = math::digamma(x.second);
                double first = lnG + digammaAlphapBeta - digammaAlpha;
                double second = lnG1m + digammaAlphapBeta - digammaBeta;
                return DoublePair(first, second);
            },
            [](DoublePair x) {
                double trigammaAlphapBeta = math::trigamma(x.first + x.second);
                double trigammaAlpha = math::trigamma(x.first);
                double trigammaBeta = math::trigamma(x.second);
                DoublePair first = std::make_pair(trigammaAlphapBeta - trigammaAlpha, trigammaAlphapBeta);
                DoublePair second = std::make_pair(trigammaAlphapBeta, trigammaAlphapBeta - trigammaBeta);
                return std::make_tuple(first, second);
            },
            shapes))
        throw this->fitError(FitErrorCode::Undefined, "Error in root-finding procedure.");
    SetShapes(shapes.first, shapes.second);
}

template <typename RealType>
void BetaRand<RealType>::FitShapes(std::span<const RealType> sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b));

    long double lnG = this->GetSampleLogMeanNorm(sample);
    if (!std::isfinite(lnG))
        throw this->fitError(FitErrorCode::WrongReturn, this->ALPHA_ZERO);
    long double lnG1m = this->GetSampleLog1mMeanNorm(sample);
    if (!std::isfinite(lnG1m))
        throw this->fitError(FitErrorCode::WrongReturn, this->BETA_ZERO);

    /// get initial values for shapes by method of moments
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    FitShapes(lnG, lnG1m, stats.first, stats.second);
}
// ARCSINERAND

template <typename RealType>
std::string ArcsineRand<RealType>::Name() const
{
    return "Arcsine(" + this->toStringWithPrecision(GetShape()) + ", " + this->toStringWithPrecision(this->MinValue()) +
           ", " + this->toStringWithPrecision(this->MaxValue()) + ")";
}

template <typename RealType>
void ArcsineRand<RealType>::SetShape(double shape)
{
    BetaDistribution<RealType>::SetShapes(1.0 - shape, shape);
}

template <typename RealType>
void ArcsineRand<RealType>::FitShape(long double lnG, long double lnG1m)
{
    double shape = randlib::consts::pi / (lnG1m - lnG);
    if (!std::isfinite(shape))
        SetShape(0.5);
    shape = -randlib::consts::inv_pi * math::atan(shape);
    SetShape(shape > 0 ? shape : shape + 1);
}

template <typename RealType>
void ArcsineRand<RealType>::FitShape(std::span<const RealType> sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b));

    int n = sample.size();
    long double lnG = 0, lnG1m = 0;
    for (double var : sample) {
        double x = (var - this->a) * this->bmaInv;
        lnG += std::log(x);
        lnG1m += std::log1pl(-x);
    }
    if (!std::isfinite(lnG))
        throw this->fitError(FitErrorCode::WrongReturn, this->ALPHA_ZERO);
    if (!std::isfinite(lnG1m))
        throw this->fitError(FitErrorCode::WrongReturn, this->BETA_ZERO);
    lnG /= n;
    lnG1m /= n;
    FitShape(lnG, lnG1m);
}
// BALDINGNICHOLSRAND

template <typename RealType>
BaldingNicholsRand<RealType>::BaldingNicholsRand(double fixatingIndex, double frequency)
{
    SetFixatingIndexAndFrequency(fixatingIndex, frequency);
}

template <typename RealType>
std::string BaldingNicholsRand<RealType>::Name() const
{
    return "Balding-Nichols(" + this->toStringWithPrecision(GetFixatingIndex()) + ", " +
           this->toStringWithPrecision(GetFrequency()) + ")";
}

template <typename RealType>
void BaldingNicholsRand<RealType>::SetFixatingIndexAndFrequency(double fixatingIndex, double frequency)
{
    if (!(fixatingIndex > 0 && fixatingIndex < 1))
        throw std::invalid_argument("Balding-Nichols distribution: fixation index should be in interval (0, 1)");
    if (!(frequency > 0 && frequency < 1))
        throw std::invalid_argument("Balding-Nichols distribution: frequency should be in interval (0, 1)");
    fixationIndex = fixatingIndex;
    p = frequency;

    double frac = (1.0 - fixationIndex) / fixationIndex, fracP = frac * p;
    BetaDistribution<RealType>::SetShapes(fracP, frac - fracP);
}

} // namespace randlib

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END

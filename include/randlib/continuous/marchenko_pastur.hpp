#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/beta.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The MarchenkoPasturRand class <BR>
 * Marchenko-Pastur distribution
 *
 * Notation: X ~ Marchenko-Pastur(λ, σ)
 */
template <typename RealType>
class MarchenkoPasturRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    double lambda = 1; ///< ratio index λ
    double sigmaSq = 1; ///< scale parameter σ^2
    double a = 0; ///< minimum value (apart from 0 if λ > 1)
    double b = 4; ///< maximum value
    double logLambda = 0; /// < log(λ)

    BetaRand<RealType> BetaRV{0.5, 1.5, 0, 4}; ///< beta-distributed rv for generator
    double M = 1.0; ///< rejection constant

public:
    MarchenkoPasturRand(double ratio = 1, double scale = 1);
    [[nodiscard]] std::string Name() const override;

    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return (lambda < 1) ? sigmaSq * a : 0; }
    [[nodiscard]] RealType MaxValue() const override { return sigmaSq * b; }

    void SetParameters(double ratio, double scale);
    [[nodiscard]] double GetRatio() const noexcept { return lambda; }
    [[nodiscard]] double GetScale() const noexcept { return sigmaSq; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;

private:
    /**
     * @param x
     * @return S(x) for λ > 1
     */
    double ccdfForLargeRatio(const RealType &x) const;
    /**
     * @param x
     * @return F(x) for 0 < λ <= 1
     */
    double cdfForSmallRatio(const RealType &x) const;

public:
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

private:
    enum GENERATOR_ID
    {
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

    template <std::uniform_random_bit_generator G>
    RealType variateForTinyRatio(G &g) const;
    template <std::uniform_random_bit_generator G>
    RealType variateForSmallRatio(G &g) const;
    template <std::uniform_random_bit_generator G>
    RealType variateForLargeRatio(G &g) const;
    template <std::uniform_random_bit_generator G>
    RealType variateForHugeRatio(G &g) const;

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

private:
    [[nodiscard]] long double Moment(int n) const;

public:
    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template <typename RealType>
MarchenkoPasturRand<RealType>::MarchenkoPasturRand(double ratio, double scale)
{
    SetParameters(ratio, scale);
}

template <typename RealType>
std::string MarchenkoPasturRand<RealType>::Name() const
{
    return "Marchenko-Pastur(" + this->toStringWithPrecision(GetRatio()) + ", " +
           this->toStringWithPrecision(GetScale()) + ")";
}

template <typename RealType>
void MarchenkoPasturRand<RealType>::SetParameters(double ratio, double scale)
{
    if (!(ratio > 0.0))
        throw std::invalid_argument("Marchenko-Pastur distribution: ratio parameter should be positive");
    if (!(scale > 0.0))
        throw std::invalid_argument("Marchenko-Pastur distribution: scale should be positive");
    lambda = ratio;
    double sqrtLambda = std::sqrt(lambda);
    a = 1.0 - sqrtLambda;
    a *= a;
    b = 1.0 + sqrtLambda;
    b *= b;

    sigmaSq = scale;
    logLambda = std::log(lambda);

    GENERATOR_ID genId = getIdOfUsedGenerator();
    /// the continuous part has density proportional to sqrt((b - x)(x - a)) / x for every ratio; M makes the
    /// largest acceptance probability one, since a smaller M only multiplies the expected number of iterations
    if (genId == TINY_RATIO || genId == HUGE_RATIO) {
        BetaRV.SetShapes(1.5, 1.5);
        BetaRV.SetSupport(a, b);
        M = a;
    } else {
        BetaRV.SetShapes(0.5, 1.5);
        BetaRV.SetSupport(0, b);
        M = 0.25 * b / sqrtLambda;
    }
}

template <typename RealType>
double MarchenkoPasturRand<RealType>::f(const RealType &x) const
{
    if (x == 0.0)
        return (lambda < 1) ? 0.0 : INFINITY;
    double xSt = x / sigmaSq;
    if (xSt < a || xSt > b)
        return 0.0;
    double y = 0.5 * std::sqrt((b - xSt) * (xSt - a));
    y /= (lambda * xSt);
    y /= (randlib::consts::pi * sigmaSq);
    return y;
}

template <typename RealType>
double MarchenkoPasturRand<RealType>::logf(const RealType &x) const
{
    if (x == 0.0)
        return (lambda < 1) ? -INFINITY : INFINITY;
    double xSt = x / sigmaSq;
    if (xSt < a || xSt > b)
        return -INFINITY;
    double y = 0.5 * std::log((b - xSt) * (xSt - a));
    y -= randlib::consts::ln2 + randlib::consts::lnpi + logLambda + std::log(x);
    return y;
}

template <typename RealType>
double MarchenkoPasturRand<RealType>::ccdfForLargeRatio(const RealType &x) const
{
    double y1 = 1.0 - x + lambda;
    double lambdam1 = lambda - 1.0;
    double lambdap1 = lambda + 1.0;
    double temp = std::sqrt(4 * lambda - y1 * y1);
    if (temp != 0.0) {
        y1 /= temp;
        y1 = lambdap1 * math::atan(y1);
    } else
        y1 = lambdap1 * math::sign(y1) * randlib::consts::half_pi;
    double y2 = x * lambdap1;
    y2 -= lambdam1 * lambdam1;
    y2 /= temp * lambdam1;
    y2 = lambdam1 * math::atan(y2);
    double y = randlib::consts::pi - temp + y1 + y2;
    y /= randlib::consts::pi * lambda;
    return 0.5 * y;
}

template <typename RealType>
double MarchenkoPasturRand<RealType>::cdfForSmallRatio(const RealType &x) const
{
    double y1 = 1.0 - x + lambda;
    double temp = std::sqrt(4 * lambda - y1 * y1);
    double lambdam1 = lambda - 1.0;
    double lambdap1 = lambda + 1.0;
    if (temp != 0.0) {
        y1 /= temp;
        y1 = lambdap1 * math::atan(y1);
    } else
        y1 = lambdap1 * math::sign(y1) * randlib::consts::half_pi;
    double y2 = 0.0;
    if (lambdam1 != 0) {
        y2 = x * lambdap1;
        y2 -= lambdam1 * lambdam1;
        y2 /= -temp * lambdam1;
        y2 = lambdam1 * math::atan(y2);
    }
    double y = randlib::consts::pi * lambda + temp - y1 + y2;
    y /= randlib::consts::pi * lambda;
    return 0.5 * y;
}

template <typename RealType>
double MarchenkoPasturRand<RealType>::F(const RealType &x) const
{
    double xSt = x / sigmaSq;
    if (xSt < 0.0)
        return 0.0;
    if (xSt >= b)
        return 1.0;
    if (lambda > 1.0)
        return (xSt > a) ? 1.0 - ccdfForLargeRatio(xSt) : 1.0 - 1.0 / lambda;
    return (xSt > a) ? cdfForSmallRatio(xSt) : 0.0;
}

template <typename RealType>
double MarchenkoPasturRand<RealType>::S(const RealType &x) const
{
    double xSt = x / sigmaSq;
    if (xSt < 0.0)
        return 1.0;
    if (xSt >= b)
        return 0.0;
    if (lambda > 1.0)
        return (xSt > a) ? ccdfForLargeRatio(xSt) : 1.0 / lambda;
    return (xSt > a) ? 1.0 - cdfForSmallRatio(xSt) : 1.0;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType MarchenkoPasturRand<RealType>::variateForTinyRatio(G &g) const
{
    size_t iter = 0;
    do {
        double X = BetaRV(g);
        double U = detail::uniform01<RealType>(g);
        if (U < M / X)
            return X;
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Marchenko-Pastur distribution: sampling failed");
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType MarchenkoPasturRand<RealType>::variateForSmallRatio(G &g) const
{
    size_t iter = 0;
    do {
        double X = BetaRV(g);
        double U = detail::uniform01<RealType>(g);
        double ratio = M * (1.0 - a / X);
        if (U * U < ratio)
            return X;
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Marchenko-Pastur distribution: sampling failed");
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType MarchenkoPasturRand<RealType>::variateForLargeRatio(G &g) const
{
    return (detail::uniform01<RealType>(g) > 1.0 / lambda) ? 0.0 : variateForSmallRatio(g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType MarchenkoPasturRand<RealType>::variateForHugeRatio(G &g) const
{
    return (detail::uniform01<RealType>(g) > 1.0 / lambda) ? 0.0 : variateForTinyRatio(g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType MarchenkoPasturRand<RealType>::operator()(G &g) const
{
    switch (getIdOfUsedGenerator()) {
    case TINY_RATIO:
        return sigmaSq * variateForTinyRatio(g);
    case SMALL_RATIO:
        return sigmaSq * variateForSmallRatio(g);
    case LARGE_RATIO:
        return sigmaSq * variateForLargeRatio(g);
    case HUGE_RATIO:
        return sigmaSq * variateForHugeRatio(g);
    default:
        throw std::runtime_error("Marchenko-Pastur distribution: invalid generator id");
    }
}

template <typename RealType>
void MarchenkoPasturRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    switch (getIdOfUsedGenerator()) {
    case TINY_RATIO:
        for (RealType &var : outputData)
            var = sigmaSq * variateForTinyRatio(g);
        break;
    case SMALL_RATIO:
        for (RealType &var : outputData)
            var = sigmaSq * variateForSmallRatio(g);
        break;
    case LARGE_RATIO:
        for (RealType &var : outputData)
            var = sigmaSq * variateForLargeRatio(g);
        break;
    case HUGE_RATIO:
        for (RealType &var : outputData)
            var = sigmaSq * variateForHugeRatio(g);
        break;
    default:
        return;
    }
}

template <typename RealType>
long double MarchenkoPasturRand<RealType>::Moment(int n) const
{
    if (!(n >= 0))
        throw std::invalid_argument("Marchenko-Pastur distribution: degree n should be non-negative");
    switch (n) {
    case 0:
        return 1.0;
    case 1:
        return sigmaSq;
    case 2:
        return sigmaSq * sigmaSq * lambda;
    default: {
        long double sum = 0.0;
        long double temp = math::lfact(n) + math::lfact(n - 1);
        long double nlogSigmaSq = n * std::log(sigmaSq);
        for (int k = 0; k != n; ++k) {
            long double term = temp;
            int kp1 = k + 1;
            term -= 2 * math::lfact(k);
            term -= math::lfact(n - k);
            term -= math::lfact(n - kp1);
            term += k * logLambda;
            term += nlogSigmaSq;
            term = std::exp(term) / kp1;
            sum += term;
        }
        return sum;
    }
    }
}

template <typename RealType>
long double MarchenkoPasturRand<RealType>::Mean() const
{
    return sigmaSq;
}

template <typename RealType>
long double MarchenkoPasturRand<RealType>::Variance() const
{
    return sigmaSq * sigmaSq * lambda;
}

template <typename RealType>
RealType MarchenkoPasturRand<RealType>::Mode() const
{
    if (lambda > 1)
        return 0.0;
    RealType mode = lambda - 1.0;
    mode *= mode;
    mode /= lambda + 1.0;
    return sigmaSq * mode;
}

template <typename RealType>
long double MarchenkoPasturRand<RealType>::Skewness() const
{
    long double mu = Mean();
    long double var = Variance();
    long double skewness = Moment(3);
    skewness -= std::pow(mu, 3);
    skewness -= 3 * mu * var;
    skewness /= std::pow(var, 1.5);
    return skewness;
}

template <typename RealType>
long double MarchenkoPasturRand<RealType>::ExcessKurtosis() const
{
    long double mu = Mean();
    long double var = Variance();
    long double moment3 = Moment(3);
    long double kurtosis = Moment(4);
    long double muSq = mu * mu;
    kurtosis -= 4 * mu * moment3;
    kurtosis += 3 * muSq * muSq;
    kurtosis += 6 * muSq * var;
    kurtosis /= var * var;
    return kurtosis - 3.0;
}

template <typename RealType>
RealType MarchenkoPasturRand<RealType>::quantileImpl(double p) const
{
    if (p <= 1.0 - 1.0 / lambda)
        return 0.0;
    RealType minInitValue = sigmaSq * a;
    RealType maxInitValue = sigmaSq * b;
    RealType initValue = minInitValue + p * (maxInitValue - minInitValue);
    return ContinuousDistribution<RealType>::quantileImpl(p, initValue);
}

template <typename RealType>
RealType MarchenkoPasturRand<RealType>::quantileImpl1m(double p) const
{
    if (p >= 1.0 / lambda)
        return 0.0;
    RealType minInitValue = sigmaSq * a;
    RealType maxInitValue = sigmaSq * b;
    RealType initValue = maxInitValue - p * (maxInitValue - minInitValue);
    return ContinuousDistribution<RealType>::quantileImpl1m(p, initValue);
}

template <typename RealType>
std::complex<double> MarchenkoPasturRand<RealType>::CFImpl(double t) const
{
    if (lambda < 1)
        return ContinuousDistribution<RealType>::CFImpl(t);
    /// otherwise we have singularity at point 0
    if (lambda == 1) {
        /// we split integrand for real part on (cos(tx)-1)f(x) and f(x)
        double re = this->ExpectedValue([t](double x) { return std::cos(t * x) - 1.0; }, 0, 4 * sigmaSq);
        double im = this->ExpectedValue([t](double x) { return std::sin(t * x); }, 0, 4 * sigmaSq);
        return std::complex<double>(1.0 + re, im);
    }
    /// for λ > 1 we split integral on 2 parts: at point 0 and the rest
    double re = this->ExpectedValue([t](double x) { return std::cos(t * x); }, sigmaSq * a, sigmaSq * b);
    double im = this->ExpectedValue([t](double x) { return std::sin(t * x); }, sigmaSq * a, sigmaSq * b);
    return std::complex<double>(1.0 - 1.0 / lambda + re, im);
}

} // namespace randlib

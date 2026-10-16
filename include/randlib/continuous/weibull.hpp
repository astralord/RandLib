#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The WeibullRand class <BR>
 * Weibull distribution
 *
 * Notation: X ~ Weibull(λ, k)
 */
template <typename RealType>
class WeibullRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

    double lambda = 1; ///< scale λ
    double k = 1; ///< shape k
    double kInv = 1; /// 1 /k
    double logk_lambda = 0; /// log(k/λ)

public:
    WeibullRand(double scale = 1, double shape = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double scale, double shape);
    [[nodiscard]] inline double GetScale() const noexcept { return lambda; }
    [[nodiscard]] inline double GetShape() const noexcept { return k; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Entropy() const;

private:
    /**
     * @brief getNorm
     * @param sample
     * @return L^k norm of the sample
     */
    double getNorm(std::span<const RealType> sample) const;

public:
    /**
     * @fn FitScale
     * Fit λ by maximum-likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitScale(std::span<const RealType> sample);

    /**
     * @fn FitScaleBayes
     * Fit λ, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER InverseGammaRand<RealType>
    FitScaleBayes(std::span<const RealType> sample, const InverseGammaRand<RealType> &priorDistribution,
                  bool MAP = false);
};

} // namespace randlib

#include <randlib/continuous/inverse_gamma.hpp>
#include <randlib/continuous/exponential.hpp>

namespace randlib
{

template <typename RealType>
WeibullRand<RealType>::WeibullRand(double scale, double shape)
{
    SetParameters(scale, shape);
}

template <typename RealType>
std::string WeibullRand<RealType>::Name() const
{
    return "Weibull(" + this->toStringWithPrecision(GetScale()) + ", " + this->toStringWithPrecision(GetShape()) + ")";
}

template <typename RealType>
void WeibullRand<RealType>::SetParameters(double scale, double shape)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Weibull distribution: scale should be positive");
    if (!(shape > 0.0))
        throw std::invalid_argument("Weibull distribution: shape should be positive");
    lambda = scale;
    k = shape;
    kInv = 1.0 / k;
    logk_lambda = std::log(k / lambda);
}

template <typename RealType>
double WeibullRand<RealType>::f(const RealType &x) const
{
    if (x < 0)
        return 0;
    if (x == 0) {
        if (k == 1)
            return 1.0 / lambda;
        return (k > 1) ? 0.0 : INFINITY;
    }
    return std::exp(this->logf(x));
}

template <typename RealType>
double WeibullRand<RealType>::logf(const RealType &x) const
{
    if (x < 0)
        return -INFINITY;
    if (x == 0) {
        if (k == 1)
            return logk_lambda;
        return (k > 1) ? -INFINITY : INFINITY;
    }
    double xAdj = x / lambda;
    double xAdjPow = std::pow(xAdj, k);
    return logk_lambda + (k - 1) * std::log(xAdj) - xAdjPow;
}

template <typename RealType>
double WeibullRand<RealType>::F(const RealType &x) const
{
    return (x > 0.0) ? -std::expm1l(-std::pow(x / lambda, k)) : 0.0;
}

template <typename RealType>
double WeibullRand<RealType>::S(const RealType &x) const
{
    return (x > 0.0) ? std::exp(-std::pow(x / lambda, k)) : 1.0;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType WeibullRand<RealType>::operator()(G &g) const
{
    return lambda * std::pow(detail::standardExponential<RealType>(g), kInv);
}

template <typename RealType>
long double WeibullRand<RealType>::Mean() const
{
    return lambda * std::tgammal(1 + kInv);
}

template <typename RealType>
long double WeibullRand<RealType>::Variance() const
{
    double res = std::tgammal(1 + kInv);
    res *= -res;
    res += std::tgammal(1 + kInv + kInv);
    return lambda * lambda * res;
}

template <typename RealType>
RealType WeibullRand<RealType>::Median() const
{
    return lambda * std::pow(randlib::consts::ln2, kInv);
}

template <typename RealType>
RealType WeibullRand<RealType>::Mode() const
{
    if (k <= 1)
        return 0;
    double y = std::log1pl(-kInv);
    y = std::exp(kInv * y);
    return lambda * y;
}

template <typename RealType>
long double WeibullRand<RealType>::Skewness() const
{
    long double mu = Mean();
    long double var = Variance();
    long double sigma = std::sqrt(var);
    long double numerator = std::tgammal(1 + 3.0 * kInv);
    numerator *= lambda * lambda * lambda;
    numerator -= 3 * mu * var;
    numerator -= mu * mu * mu;
    long double denominator = var * sigma;
    return numerator / denominator;
}

template <typename RealType>
long double WeibullRand<RealType>::ExcessKurtosis() const
{
    long double mu = Mean();
    long double var = Variance();
    long double sigma = std::sqrt(var);
    long double skewness = Skewness();
    long double numerator = lambda * lambda;
    numerator *= numerator;
    numerator *= std::tgammal(1 + 4.0 * kInv);
    numerator -= 4 * skewness * var * sigma * mu;
    long double mu2 = mu * mu;
    numerator -= 6 * mu2 * var;
    numerator -= mu2 * mu2;
    long double kurtosis = numerator / (var * var);
    return kurtosis - 3;
}

template <typename RealType>
RealType WeibullRand<RealType>::quantileImpl(double p) const
{
    double x = -std::log1pl(-p);
    x = std::pow(x, kInv);
    return lambda * x;
}

template <typename RealType>
RealType WeibullRand<RealType>::quantileImpl1m(double p) const
{
    double x = -std::log(p);
    x = std::pow(x, kInv);
    return lambda * x;
}

template <typename RealType>
std::complex<double> WeibullRand<RealType>::CFImpl(double t) const
{
    double lambdaT = lambda * t;
    if (k >= 1) {
        if (lambdaT > 0.5)
            return ContinuousDistribution<RealType>::CFImpl(t);
        /// for λt < 0.5, the worst case scenario for series expansion is n ~ 70
        long double re = 0.0, im = 0.0;
        long double addon = 0.0;
        long double logLambdaT = std::log(lambdaT);
        /// Series representation for real part
        int n = 0;
        do {
            int n2 = n + n;
            addon = n2 * logLambdaT;
            addon += std::lgammal(1.0 + n2 / k);
            addon -= std::lgammal(1.0 + n2);
            addon = std::exp(addon);
            re += (n & 1) ? -addon : addon;
            ++n;
        } while (std::fabs(addon) > math::detail::minimumTolerance * std::fabs(re));
        /// Series representation for imaginary part
        n = 0;
        do {
            int n2p1 = n + n + 1;
            addon = n2p1 * logLambdaT;
            addon += std::lgammal(1.0 + n2p1 / k);
            addon -= std::lgammal(1.0 + n2p1);
            addon = std::exp(addon);
            im += (n & 1) ? -addon : addon;
            ++n;
        } while (std::fabs(addon) > math::detail::minimumTolerance * std::fabs(im));
        return std::complex<double>(re, im);
    }

    /// For real part with k < 1 we split the integral on two intervals
    double re1 = math::integrate(
        [this, t](double x) {
            if (x <= 0.0 || x > 1.0)
                return 0.0;
            double xAdj = x / lambda;
            double xAdjPow = std::pow(xAdj, k - 1);
            double y = k / lambda * xAdjPow * std::expm1l(-xAdj * xAdjPow);
            return std::cos(t * x) * y;
        },
        0.0, 1.0);

    double re2 = this->ExpectedValue([t](double x) { return std::cos(t * x); }, 1.0, INFINITY);

    double re3 = t * math::integrate(
                         [this, t](double x) {
                             if (x <= 0.0)
                                 return 0.0;
                             return std::sin(t * x) * std::pow(x, k);
                         },
                         0.0, 1.0);
    re3 += std::cos(t);
    re3 /= std::pow(lambda, k);

    double re = re1 + re2 + re3;

    double im = this->ExpectedValue([t](double x) { return std::sin(t * x); }, 0.0, INFINITY);

    return std::complex<double>(re, im);
}

template <typename RealType>
long double WeibullRand<RealType>::Entropy() const
{
    return randlib::consts::egamma * (1.0 - kInv) + std::log(lambda * kInv) + 1.0;
}

template <typename RealType>
double WeibullRand<RealType>::getNorm(std::span<const RealType> sample) const
{
    long double sum = 0;
    RealType maxVar = *std::max_element(sample.begin(), sample.end());
    for (const RealType &var : sample) {
        sum += std::pow(var / maxVar, k);
    }
    long double avg = sum / sample.size();
    return maxVar * std::pow(avg, kInv);
}

template <typename RealType>
void WeibullRand<RealType>::FitScale(std::span<const RealType> sample)
{
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    double knorm = getNorm(sample);
    SetParameters(knorm, k);
}

template <typename RealType>
InverseGammaRand<RealType> WeibullRand<RealType>::FitScaleBayes(std::span<const RealType> sample,
                                                                const InverseGammaRand<RealType> &priorDistribution,
                                                                bool MAP)
{
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    int n = sample.size();
    double newShape = priorDistribution.GetShape() + n;
    double knorm = getNorm(sample);
    double newRate = priorDistribution.GetRate() + n * std::pow(knorm, k);
    InverseGammaRand<RealType> posteriorDistribution(newShape, newRate);
    double powScale = MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean();
    SetParameters(std::pow(powScale, kInv), k);
    return posteriorDistribution;
}

} // namespace randlib

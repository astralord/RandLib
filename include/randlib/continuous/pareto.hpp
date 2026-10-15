#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The ParetoRand class <BR>
 * Pareto distribution
 *
 * Notation: X ~ Pareto(α, σ)
 *
 * Related distributions: <BR>
 * ln(X / σ) ~ Exp(α)
 */
template <typename RealType>
class ParetoRand final : public ContinuousDistribution<RealType>
{
    double alpha = 1; ///< shape α
    double sigma = 1; ///< scale σ
    double logAlpha = 0; ///< log(α)
    double logSigma = 0; ///< log(σ)

public:
    ParetoRand(double shape = 1, double scale = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return sigma; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetShape(double shape);
    void SetScale(double scale);
    [[nodiscard]] inline double GetShape() const noexcept { return alpha; }
    [[nodiscard]] inline double GetScale() const noexcept { return sigma; }
    [[nodiscard]] inline double GetLogShape() const noexcept { return logAlpha; }
    [[nodiscard]] inline double GetLogScale() const noexcept { return logSigma; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

private:
    template <std::uniform_random_bit_generator G>
    static RealType variateForAlphaEqualOne(G &g);
    template <std::uniform_random_bit_generator G>
    static RealType variateForAlphaEqualTwo(G &g);
    template <std::uniform_random_bit_generator G>
    static RealType variateForGeneralAlpha(double shape, G &g);

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(double shape, G &g);

protected:
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

public:
    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

public:
    [[nodiscard]] long double Entropy() const;

    /**
     * @fn FitShape
     * Fit α to maximum-likelihood estimator
     * or to UMVU if unbiased == true
     * @param sample
     * @param unbiased
     */
    void FitShape(std::span<const RealType> sample, bool unbiased = false);

    /**
     * @fn FitScale
     * Fit σ to maximum-likelihood estimator
     * or to UMVU if unbiased == true
     * @param sample
     * @param unbiased
     */
    void FitScale(std::span<const RealType> sample, bool unbiased = false);

    /**
     * @fn Fit
     * Fit parameters to maximum-likelihood estimators
     * or to UMVU if unbiased == true
     * @param sample
     * @param unbiased
     */
    void Fit(std::span<const RealType> sample, bool unbiased = false);

    /**
     * @fn FitShapeBayes
     * Fit α, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution of α
     */
    GammaRand<RealType> FitShapeBayes(std::span<const RealType> sample,
                                      const GammaDistribution<RealType> &priorDistribution, bool MAP = false);
};

} // namespace randlib

#include <randlib/continuous/exponential.hpp>
#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template <typename RealType>
ParetoRand<RealType>::ParetoRand(double shape, double scale)
{
    SetShape(shape);
    SetScale(scale);
}

template <typename RealType>
std::string ParetoRand<RealType>::Name() const
{
    return "Pareto(" + this->toStringWithPrecision(GetShape()) + ", " + this->toStringWithPrecision(GetScale()) + ")";
}

template <typename RealType>
void ParetoRand<RealType>::SetShape(double shape)
{
    if (!(shape > 0.0))
        throw std::invalid_argument("Pareto distribution: shape should be positive");
    alpha = shape;
    logAlpha = std::log(alpha);
}

template <typename RealType>
void ParetoRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Pareto distribution: scale should be positive");
    sigma = scale;
    logSigma = std::log(sigma);
}

template <typename RealType>
double ParetoRand<RealType>::f(const RealType &x) const
{
    return (x < sigma) ? 0.0 : std::exp(logf(x));
}

template <typename RealType>
double ParetoRand<RealType>::logf(const RealType &x) const
{
    if (x < sigma)
        return -INFINITY;
    double logX = std::log(x);
    double y = logSigma - logX;
    y *= alpha;
    y -= logX;
    y += logAlpha;
    return y;
}

template <typename RealType>
double ParetoRand<RealType>::F(const RealType &x) const
{
    return (x > sigma) ? -std::expm1l(alpha * std::log(sigma / x)) : 0.0;
}

template <typename RealType>
double ParetoRand<RealType>::S(const RealType &x) const
{
    return (x > sigma) ? std::pow(sigma / x, alpha) : 1.0;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType ParetoRand<RealType>::variateForAlphaEqualOne(G &g)
{
    return 1.0 / detail::uniform01<RealType>(g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType ParetoRand<RealType>::variateForAlphaEqualTwo(G &g)
{
    return 1.0 / std::sqrt(detail::uniform01<RealType>(g));
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType ParetoRand<RealType>::variateForGeneralAlpha(double shape, G &g)
{
    return std::exp(detail::standardExponential<RealType>(g) / shape);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType ParetoRand<RealType>::operator()(G &g) const
{
    return sigma * StandardVariate(alpha, g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType ParetoRand<RealType>::StandardVariate(double shape, G &g)
{
    if (math::areClose(shape, 1.0))
        return variateForAlphaEqualOne(g);
    if (math::areClose(shape, 2.0))
        return variateForAlphaEqualTwo(g);
    return variateForGeneralAlpha(shape, g);
}

template <typename RealType>
void ParetoRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    if (math::areClose(alpha, 1.0)) {
        for (RealType &var : outputData)
            var = sigma * variateForAlphaEqualOne(g);
    } else if (math::areClose(alpha, 2.0)) {
        for (RealType &var : outputData)
            var = sigma * variateForAlphaEqualTwo(g);
    } else {
        for (RealType &var : outputData)
            var = sigma * variateForGeneralAlpha(alpha, g);
    }
}

template <typename RealType>
long double ParetoRand<RealType>::Mean() const
{
    return (alpha > 1) ? alpha * sigma / (alpha - 1) : INFINITY;
}

template <typename RealType>
long double ParetoRand<RealType>::Variance() const
{
    if (alpha > 2) {
        long double var = sigma / (alpha - 1);
        var *= var;
        return alpha * var / (alpha - 2);
    }
    return (alpha > 1) ? INFINITY : NAN;
}

template <typename RealType>
RealType ParetoRand<RealType>::Median() const
{
    return std::exp(logSigma + randlib::consts::ln2 / alpha);
}

template <typename RealType>
RealType ParetoRand<RealType>::Mode() const
{
    return sigma;
}

template <typename RealType>
long double ParetoRand<RealType>::Skewness() const
{
    if (alpha <= 3)
        return INFINITY;
    double skewness = (alpha - 2.0) / alpha;
    skewness = std::sqrt(skewness);
    skewness *= (1 + alpha) / (alpha - 3);
    return skewness + skewness;
}

template <typename RealType>
long double ParetoRand<RealType>::ExcessKurtosis() const
{
    if (alpha <= 4)
        return INFINITY;
    double numerator = alpha + 1;
    numerator *= alpha;
    numerator -= 6;
    numerator *= alpha;
    numerator -= 2;
    double denominator = alpha * (alpha - 3) * (alpha - 4);
    return 6.0 * numerator / denominator;
}

template <typename RealType>
RealType ParetoRand<RealType>::quantileImpl(double p) const
{
    double y = logSigma - std::log1pl(-p) / alpha;
    return std::exp(y);
}

template <typename RealType>
RealType ParetoRand<RealType>::quantileImpl1m(double p) const
{
    double y = logSigma - std::log(p) / alpha;
    return std::exp(y);
}

template <typename RealType>
long double ParetoRand<RealType>::Entropy() const
{
    return logSigma - logAlpha + 1.0 / alpha + 1;
}

template <typename RealType>
void ParetoRand<RealType>::FitShape(std::span<const RealType> sample, bool unbiased)
{
    if (!this->allElementsAreNotSmallerThan(sigma, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(sigma));
    double invShape = this->GetSampleLogMean(sample) - logSigma;
    if (invShape == 0.0)
        throw this->fitError(FitErrorCode::WrongSample,
                             "Possibly all the elements of the sample coincide with the lower boundary σ.");
    double shape = 1.0 / invShape;
    if (unbiased)
        shape *= (1.0 - 1.0 / sample.size());
    SetShape(shape);
}

template <typename RealType>
void ParetoRand<RealType>::FitScale(std::span<const RealType> sample, bool unbiased)
{
    RealType minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    double scale = unbiased ? (1.0 - 1.0 / (sample.size() * alpha)) * minVar : minVar;
    SetScale(scale);
}

template <typename RealType>
void ParetoRand<RealType>::Fit(std::span<const RealType> sample, bool unbiased)
{
    RealType minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);

    double logBiasedSigma = std::log(minVar);
    double invShape = this->GetSampleLogMean(sample) - logBiasedSigma;
    if (invShape == 0.0)
        throw this->fitError(FitErrorCode::WrongSample, "Possibly all the elements of the sample are the same.");
    double shape = 1.0 / invShape;

    double scale = minVar;
    if (unbiased) {
        int n = sample.size();
        shape *= 1.0 - 2.0 / n;
        scale *= 1.0 - invShape / (n - 1);
    }
    SetScale(scale);
    SetShape(shape);
}

template <typename RealType>
GammaRand<RealType> ParetoRand<RealType>::FitShapeBayes(std::span<const RealType> sample,
                                                        const GammaDistribution<RealType> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNotSmallerThan(sigma, sample))
        throw this->fitError(FitErrorCode::WrongSample,
                             this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(sigma));
    int n = sample.size();
    double newShape = priorDistribution.GetShape() + n;
    double newRate = priorDistribution.GetRate() + n * (this->GetSampleLogMean(sample) - logSigma);
    GammaRand<RealType> posteriorDistribution(newShape, newRate);
    SetShape(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

} // namespace randlib

#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/beta.hpp>
#include <randlib/core/continuous.hpp>

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN

namespace randlib
{

/**
 * @brief The BetaPrimeRand class <BR>
 * Beta-prime distribution
 *
 * f(x | α, β) = x^{α-1} (1 + x)^{-α - β} / B(α, β), <BR>
 * where B(α, β) denotes Beta function
 *
 * Notation: X ~ B'(α, β)
 *
 * Related distributions: <BR>
 * X / (X + 1) ~ B(α, β) <BR>
 * X = Y / Z, where Y ~ Γ(α) and Z ~ Γ(β) <BR>
 * β/α * X ~ F(2α, 2β)
 */
template <typename RealType>
class BetaPrimeRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    double alpha = 1; ///< first shape α
    double beta = 1; ///< second shape β
    BetaRand<RealType> B{};

public:
    BetaPrimeRand(double shape1 = 1, double shape2 = 1);
    [[nodiscard]] std::string Name() const override;
    void SetShapes(double shape1, double shape2);
    [[nodiscard]] inline double GetAlpha() const noexcept { return alpha; }
    [[nodiscard]] inline double GetBeta() const noexcept { return beta; }

    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

private:
    RealType fromBetaVariate(const RealType &betaVar) const;

public:
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
    /**
     * @return B(α, β)
     */
    [[nodiscard]] inline double GetBetaFunction() const { return B.GetBetaFunction(); }

    /**
     * @return log(B(α, β))
     */
    [[nodiscard]] inline double GetLogBetaFunction() const { return B.GetLogBetaFunction(); }

    /**
     * fit α by maximum-likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitAlpha(std::span<const RealType> sample);

    /**
     * fit β by maximum-likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitBeta(std::span<const RealType> sample);

    /**
     * fit shapes by maximum-likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void Fit(std::span<const RealType> sample);
};

} // namespace randlib

namespace randlib
{

template <typename RealType>
BetaPrimeRand<RealType>::BetaPrimeRand(double shape1, double shape2)
{
    SetShapes(shape1, shape2);
}

template <typename RealType>
std::string BetaPrimeRand<RealType>::Name() const
{
    return "Beta Prime(" + this->toStringWithPrecision(GetAlpha()) + ", " + this->toStringWithPrecision(GetBeta()) +
           ")";
}

template <typename RealType>
void BetaPrimeRand<RealType>::SetShapes(double shape1, double shape2)
{
    if (!(shape1 > 0 && shape2 > 0))
        throw std::invalid_argument("Beta-prime distribution: shapes should be positive");
    B.SetShapes(shape1, shape2);
    alpha = B.GetAlpha();
    beta = B.GetBeta();
}

template <typename RealType>
double BetaPrimeRand<RealType>::f(const RealType &x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (alpha == 1.0)
            return 1.0 / GetBetaFunction();
        return (alpha > 1) ? 0.0 : INFINITY;
    }
    return std::exp(logf(x));
}

template <typename RealType>
double BetaPrimeRand<RealType>::logf(const RealType &x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0) {
        if (alpha == 1.0)
            return -GetLogBetaFunction();
        return (alpha > 1) ? -INFINITY : INFINITY;
    }
    double y = (alpha - 1) * std::log(x);
    y -= (alpha + beta) * std::log1pl(x);
    return y - GetLogBetaFunction();
}

template <typename RealType>
double BetaPrimeRand<RealType>::F(const RealType &x) const
{
    return (x > 0) ? B.F(x / (1.0 + x)) : 0;
}

template <typename RealType>
double BetaPrimeRand<RealType>::S(const RealType &x) const
{
    if (x <= 0)
        return 1;
    /// I_{1/(1+x)}(β, α) directly: forming x / (1 + x) would round to 1 for large x
    double y = 1.0 / (1.0 + x);
    return math::ibeta(y, beta, alpha, GetLogBetaFunction(), -std::log1pl(x), std::log(x) - std::log1pl(x));
}

template <typename RealType>
RealType BetaPrimeRand<RealType>::fromBetaVariate(const RealType &betaVar) const
{
    if (betaVar > 1e-5)
        return betaVar / (1.0 - betaVar);
    RealType logVar = std::log(betaVar), log1mVar = std::log1p(-betaVar);
    return std::exp(logVar - log1mVar);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType BetaPrimeRand<RealType>::operator()(G &g) const
{
    double x = B(g);
    return fromBetaVariate(x);
}

template <typename RealType>
void BetaPrimeRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    B.Sample(g, outputData);
    for (RealType &var : outputData)
        var = fromBetaVariate(var);
}

template <typename RealType>
long double BetaPrimeRand<RealType>::Mean() const
{
    return (beta > 1) ? alpha / (beta - 1) : INFINITY;
}

template <typename RealType>
long double BetaPrimeRand<RealType>::Variance() const
{
    if (beta <= 2)
        return INFINITY;
    double betam1 = beta - 1;
    double numerator = alpha * (alpha + betam1);
    double denominator = (betam1 - 1) * betam1 * betam1;
    return numerator / denominator;
}

template <typename RealType>
RealType BetaPrimeRand<RealType>::Median() const
{
    return (alpha == beta) ? 1.0 : quantileImpl(0.5);
}

template <typename RealType>
RealType BetaPrimeRand<RealType>::Mode() const
{
    return (alpha < 1) ? 0 : (alpha - 1) / (beta + 1);
}

template <typename RealType>
long double BetaPrimeRand<RealType>::Skewness() const
{
    if (beta <= 3)
        return INFINITY;
    long double aux = alpha + beta - 1;
    long double skewness = (beta - 2) / (alpha * aux);
    skewness = std::sqrt(skewness);
    aux += alpha;
    aux += aux;
    return aux * skewness / (beta - 3);
}

template <typename RealType>
long double BetaPrimeRand<RealType>::ExcessKurtosis() const
{
    if (beta <= 4)
        return INFINITY;
    long double betam1 = beta - 1;
    long double numerator = betam1 * betam1 * (beta - 2) / (alpha * (alpha + betam1));
    numerator += 5 * beta - 11;
    long double denominator = (beta - 3) * (beta - 4);
    return 6 * numerator / denominator;
}

template <typename RealType>
RealType BetaPrimeRand<RealType>::quantileImpl(double p) const
{
    /// y / (1 - y) for the beta quantile y is only a guess: it loses accuracy as y → 1
    double y = B.Quantile(p);
    return ContinuousDistribution<RealType>::quantileImpl(p, y / (1.0 - y));
}

template <typename RealType>
RealType BetaPrimeRand<RealType>::quantileImpl1m(double p) const
{
    double y = B.Quantile1m(p);
    return ContinuousDistribution<RealType>::quantileImpl1m(p, y / (1.0 - y));
}

template <typename RealType>
std::complex<double> BetaPrimeRand<RealType>::CFImpl(double t) const
{
    return this->integrateCharacteristicFunction(t);
}

template <typename RealType>
void BetaPrimeRand<RealType>::FitAlpha(std::span<const RealType> sample)
{
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    long double lnG = this->GetSampleLogMean(sample) - B.GetSampleLog1pMeanNorm(sample);
    long double mean = 0.5;
    if (beta != 1.0) {
        mean = 0.0;
        for (const RealType &var : sample)
            mean += var / (1.0 + var);
        mean /= sample.size();
    }
    B.FitAlpha(lnG, mean);
    SetShapes(B.GetAlpha(), beta);
}

template <typename RealType>
void BetaPrimeRand<RealType>::FitBeta(std::span<const RealType> sample)
{
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    long double lnG1m = -B.GetSampleLog1pMeanNorm(sample);
    long double mean = 0.5;
    if (alpha != 1.0) {
        mean = 0.0;
        for (const RealType &var : sample)
            mean += var / (1.0 + var);
        mean /= sample.size();
    }
    B.FitBeta(lnG1m, mean);
    SetShapes(alpha, B.GetBeta());
}

template <typename RealType>
void BetaPrimeRand<RealType>::Fit(std::span<const RealType> sample)
{
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    long double lnG1m = -B.GetSampleLog1pMeanNorm(sample);
    long double lnG = this->GetSampleLogMean(sample) + lnG1m;
    long double m = 0.0, v = 0.0;
    int n = sample.size();
    for (int i = 0; i < n; ++i) {
        double x = sample[i] / (1.0 + sample[i]);
        double diff = x - m;
        m += diff / (i + 1);
        v += diff * (x - m);
    }
    B.FitShapes(lnG, lnG1m, m, v / n);
    SetShapes(B.GetAlpha(), B.GetBeta());
}

} // namespace randlib

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END

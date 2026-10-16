#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/normal.hpp>

namespace randlib
{

/**
 * @brief The LogNormalRand class <BR>
 * Log-Normal distribution
 *
 * Notation X ~ Log-Normal(μ, σ)
 *
 * Related distributions: <BR>
 * ln(X) ~ Normal(μ, σ)
 */
template <typename RealType>
class LogNormalRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

    NormalRand<RealType> X{};
    double expMu = 1; ///< exp(μ)
    double expHalfSigmaSq = 1.6487212707; ///< exp(σ^2 / 2)

public:
    LogNormalRand(double location = 0, double squaredScale = 1);

    /**
     * @fn FromLogMeanLogStd
     * @param logMean μ, mean of log(X)
     * @param logStandardDeviation σ, standard deviation of log(X)
     * @return LogNormal(μ, σ^2)
     */
    [[nodiscard]] static LogNormalRand FromLogMeanLogStd(double logMean, double logStandardDeviation)
    {
        if (!(logStandardDeviation > 0.0))
            throw std::invalid_argument("Log-Normal distribution: scale should be positive");
        return LogNormalRand(logMean, logStandardDeviation * logStandardDeviation);
    }

    /**
     * @fn FromLogMeanLogVariance
     * @param logMean μ, mean of log(X)
     * @param logVariance σ^2, variance of log(X)
     * @return LogNormal(μ, σ^2)
     */
    [[nodiscard]] static LogNormalRand FromLogMeanLogVariance(double logMean, double logVariance)
    {
        return LogNormalRand(logMean, logVariance);
    }

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    [[nodiscard]] inline double GetLocation() const { return X.Mean(); }
    [[nodiscard]] inline double GetScale() const { return X.GetScale(); }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);

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
    /**
     * @fn FitLocation
     * Fit location via maximum-likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitLocation(std::span<const RealType> sample);
    /**
     * @fn FitScale
     * Fit scale via maximum-likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitScale(std::span<const RealType> sample);
    /**
     * @fn Fit
     * Fit parameters via maximum-likelihood
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void Fit(std::span<const RealType> sample);

    /**
     * @fn FitLocationBayes
     * Set location, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER NormalRand<RealType>
    FitLocationBayes(std::span<const RealType> sample, const NormalRand<RealType> &priorDistribution, bool MAP = false);

    /**
     * @fn FitScaleBayes
     * Set scale, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER InverseGammaRand<RealType>
    FitScaleBayes(std::span<const RealType> sample, const InverseGammaRand<RealType> &priorDistribution,
                  bool MAP = false);

    /**
     * @fn FitBayes
     * Set parameters, using bayesian inference
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER NormalInverseGammaRand<RealType>
    FitBayes(std::span<const RealType> sample, const NormalInverseGammaRand<RealType> &priorDistribution,
             bool MAP = false);
};

} // namespace randlib

namespace randlib
{

template <typename RealType>
LogNormalRand<RealType>::LogNormalRand(double location, double squaredScale)
{
    if (!(squaredScale > 0.0))
        throw std::invalid_argument("Log-Normal distribution: squared scale should be positive");
    SetLocation(location);
    SetScale(std::sqrt(squaredScale));
}

template <typename RealType>
std::string LogNormalRand<RealType>::Name() const
{
    return "Log-Normal(" + this->toStringWithPrecision(GetLocation()) + ", " + this->toStringWithPrecision(GetScale()) +
           ")";
}

template <typename RealType>
void LogNormalRand<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Log-Normal distribution: location should be finite");
    X.SetLocation(location);
    expMu = std::exp(X.GetLocation());
}

template <typename RealType>
void LogNormalRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Log-Normal distribution: scale should be positive");
    X.SetScale(scale);
    expHalfSigmaSq = std::exp(0.5 * X.Variance());
}

template <typename RealType>
double LogNormalRand<RealType>::f(const RealType &x) const
{
    return (x > 0.0) ? std::exp(logf(x)) : 0.0;
}

template <typename RealType>
double LogNormalRand<RealType>::logf(const RealType &x) const
{
    if (x <= 0.0)
        return -INFINITY;
    double logX = std::log(x);
    double y = X.logf(logX);
    return y - logX;
}

template <typename RealType>
double LogNormalRand<RealType>::F(const RealType &x) const
{
    return (x > 0.0) ? X.F(std::log(x)) : 0.0;
}

template <typename RealType>
double LogNormalRand<RealType>::S(const RealType &x) const
{
    return (x > 0.0) ? X.S(std::log(x)) : 1.0;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType LogNormalRand<RealType>::operator()(G &g) const
{
    return std::exp(X(g));
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType LogNormalRand<RealType>::StandardVariate(G &g)
{
    return std::exp(detail::standardNormal<RealType>(g));
}

template <typename RealType>
long double LogNormalRand<RealType>::Mean() const
{
    return expMu * expHalfSigmaSq;
}

template <typename RealType>
long double LogNormalRand<RealType>::Variance() const
{
    double y = expMu * expHalfSigmaSq;
    return y * y * std::expm1l(X.Variance());
}

template <typename RealType>
RealType LogNormalRand<RealType>::quantileImpl(double p) const
{
    return std::exp(X.Quantile(p));
}

template <typename RealType>
RealType LogNormalRand<RealType>::quantileImpl1m(double p) const
{
    return std::exp(X.Quantile1m(p));
}

template <typename RealType>
RealType LogNormalRand<RealType>::Median() const
{
    return expMu;
}

template <typename RealType>
RealType LogNormalRand<RealType>::Mode() const
{
    return expMu / (expHalfSigmaSq * expHalfSigmaSq);
}

template <typename RealType>
long double LogNormalRand<RealType>::Skewness() const
{
    double y = std::expm1l(X.Variance());
    return (expHalfSigmaSq * expHalfSigmaSq + 2) * std::sqrt(y);
}

template <typename RealType>
long double LogNormalRand<RealType>::ExcessKurtosis() const
{
    double temp = expHalfSigmaSq * expHalfSigmaSq;
    double c = temp * temp;
    double b = c * temp;
    double a = c * c;
    return a + 2 * b + 3 * c - 6;
}

template <typename RealType>
void LogNormalRand<RealType>::FitLocation(std::span<const RealType> sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    this->SetLocation(this->GetSampleLogMean(sample));
}

template <typename RealType>
void LogNormalRand<RealType>::FitScale(std::span<const RealType> sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    double mu = X.GetLocation();
    this->SetScale(std::sqrt(this->GetSampleLogVariance(sample, mu)));
}

template <typename RealType>
void LogNormalRand<RealType>::Fit(std::span<const RealType> sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    size_t n = sample.size();
    long double logMean = 0.0L;
    long double logSqDev = 0.0L;
    for (double var : sample) {
        double logVar = std::log(var);
        logMean += logVar;
        logSqDev += logVar * logVar;
    }
    logMean /= n;
    logSqDev /= n;
    logSqDev -= logMean * logMean;

    SetLocation(logMean);
    SetScale(std::sqrt(logSqDev));
}

template <typename RealType>
NormalRand<RealType> LogNormalRand<RealType>::FitLocationBayes(std::span<const RealType> sample,
                                                               const NormalRand<RealType> &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    size_t n = sample.size();
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = X.GetPrecision();
    double numerator = n * this->GetSampleLogMean(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    NormalRand<RealType> posteriorDistribution(numerator / denominator, 1.0 / denominator);
    SetLocation(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template <typename RealType>
InverseGammaRand<RealType> LogNormalRand<RealType>::FitScaleBayes(std::span<const RealType> sample,
                                                                  const InverseGammaRand<RealType> &priorDistribution,
                                                                  bool MAP)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double newAlpha = alpha + 0.5 * n;
    double mu = X.GetLocation();
    double newBeta = beta + 0.5 * n * this->GetSampleLogVariance(sample, mu);
    InverseGammaRand<RealType> posteriorDistribution(newAlpha, newBeta);
    SetScale(std::sqrt(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean()));
    return posteriorDistribution;
}

template <typename RealType>
NormalInverseGammaRand<RealType>
LogNormalRand<RealType>::FitBayes(std::span<const RealType> sample,
                                  const NormalInverseGammaRand<RealType> &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double mu0 = priorDistribution.GetLocation();
    double lambda = priorDistribution.GetPrecision();
    DoublePair logStats = this->GetSampleLogMeanAndVariance(sample);
    double average = logStats.first, sum = n * average;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = logStats.second;
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    NormalInverseGammaRand<RealType> posteriorDistribution(newMu0, newLambda, newAlpha, newBeta);
    DoublePair newParams = MAP ? static_cast<DoublePair>(posteriorDistribution.Mode())
                               : static_cast<DoublePair>(posteriorDistribution.Mean());
    SetLocation(newParams.first);
    SetScale(std::sqrt(newParams.second));
    return posteriorDistribution;
}

} // namespace randlib

#pragma once

#include <span>

#include <array>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/exponential_family.hpp>
#include <randlib/continuous/stable.hpp>

namespace randlib
{

/**
 * @brief The NormalRand class <BR>
 * Normal distribution
 *
 * f(x | μ, σ^2) = 1 / ((2 π σ^2)^(1/2) * exp(-(x - μ)^2 / (2 σ^2))
 *
 * Notation: X ~ N(μ, σ^2)
 *
 * Related distributions: <BR>
 * X ~ S(2, 0, σ/√2, μ)
 */
template <typename RealType>
class NormalRand final : public StableDistribution<RealType>, public ExponentialFamily<RealType, DoublePair>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

    double sigma = 1; ///< scale σ

public:
    NormalRand(double location = 0, double variance = 1);

    /**
     * @fn FromMeanStd
     * @param mean μ
     * @param standardDeviation σ
     * @return N(μ, σ^2)
     */
    [[nodiscard]] static NormalRand FromMeanStd(double mean, double standardDeviation)
    {
        NormalRand X(mean);
        X.SetScale(standardDeviation);
        return X;
    }

    /**
     * @fn FromMeanVariance
     * @param mean μ
     * @param variance σ^2
     * @return N(μ, σ^2)
     */
    [[nodiscard]] static NormalRand FromMeanVariance(double mean, double variance)
    {
        return NormalRand(mean, variance);
    }
    [[nodiscard]] std::string Name() const override;

    void SetScale(double scale);
    void SetVariance(double variance);
    /**
     * @fn GetScale
     * @return σ
     */
    [[nodiscard]] inline double GetScale() const noexcept { return sigma; }
    /**
     * @fn GetLogScale
     * @return log(σ)
     */
    [[nodiscard]] inline double GetLogScale() const
    {
        return StableDistribution<RealType>::GetLogScale() - 0.5 * randlib::consts::ln2;
    }
    /**
     * @fn GetPrecision
     * @return 1/σ^2
     */
    [[nodiscard]] inline double GetPrecision() const { return 1.0 / this->Variance(); }

    DoublePair SufficientStatistic(RealType x) const override;
    DoublePair SourceParameters() const override;
    DoublePair SourceToNatural(DoublePair sourceParameters) const override;
    double LogNormalizer(DoublePair theta) const override;
    DoublePair LogNormalizerGradient(DoublePair theta) const override;
    double CarrierMeasure(RealType) const override;

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Moment(size_t n) const;
    [[nodiscard]] long double ThirdMoment() const override { return Moment(3); }
    [[nodiscard]] long double FourthMoment() const override { return Moment(4); }

    /**
     * @fn FitLocation
     * set location, returned by maximium-likelihood estimator
     * @param sample
     */
    void FitLocation(std::span<const RealType> sample);

    /**
     * @fn FitLocation
     * set location, returned by maximium-likelihood estimator
     * and return confidenceInterval for given significance level
     * @param sample
     * @param confidenceInterval
     * @param significanceLevel
     */
    void FitLocation(std::span<const RealType> sample, DoublePair &confidenceInterval, double significanceLevel);

    /**
     * @fn FitVariance
     * set variance, returned by maximium-likelihood estimator
     * @param sample
     */
    void FitVariance(std::span<const RealType> sample);

    /**
     * @fn FitVariance
     * @param sample
     * @param confidenceInterval
     * @param significanceLevel
     * @param unbiased
     */
    void FitVariance(std::span<const RealType> sample, DoublePair &confidenceInterval, double significanceLevel,
                     bool unbiased = false);

    /**
     * @fn FitScale
     * set scale, returned via maximum-likelihood estimation or unbiased estimator
     * (which might be different from the unbiased estimator of variance)
     * @param sample
     * @param unbiased
     */
    void FitScale(std::span<const RealType> sample, bool unbiased = false);

    /**
     * @fn Fit
     * set parameters, returned by maximium-likelihood estimator if unbiased = false,
     * otherwise set parameters via UMVU estimator
     * @param sample
     * @param unbiased
     */
    void Fit(std::span<const RealType> sample, bool unbiased = false);

    /**
     * @fn Fit
     * set parameters, returned by maximium-likelihood estimator if unbiased = false,
     * otherwise set parameters via UMVU estimator, and return confidence intervals for given significance level
     * @param sample
     * @param confidenceIntervalForMean
     * @param confidenceIntervalForVariance
     * @param significanceLevel
     * @param unbiased
     */
    void Fit(std::span<const RealType> sample, DoublePair &confidenceIntervalForMean,
             DoublePair &confidenceIntervalForVariance, double significanceLevel, bool unbiased = false);

    /**
     * @fn FitLocationBayes
     * set location, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    NormalRand<RealType> FitLocationBayes(std::span<const RealType> sample,
                                          const NormalRand<RealType> &priorDistribution, bool MAP = false);

    /**
     * @fn FitVarianceBayes
     * set variance, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    InverseGammaRand<RealType> FitVarianceBayes(std::span<const RealType> sample,
                                                const InverseGammaRand<RealType> &priorDistribution, bool MAP = false);

    /**
     * @fn FitBayes
     * set parameters, returned by bayesian estimation
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    NormalInverseGammaRand<RealType> FitBayes(std::span<const RealType> sample,
                                              const NormalInverseGammaRand<RealType> &priorDistribution,
                                              bool MAP = false);
};

} // namespace randlib

#include <randlib/continuous/inverse_gamma.hpp>
#include <randlib/bivariate/normal_inverse_gamma.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/student_t.hpp>

namespace randlib
{

template <typename RealType>
NormalRand<RealType>::NormalRand(double location, double variance)
    : StableDistribution<RealType>(2.0, 0.0, 1.0, location)
{
    SetVariance(variance);
}

template <typename RealType>
std::string NormalRand<RealType>::Name() const
{
    return "Normal(" + this->toStringWithPrecision(this->GetLocation()) + ", " +
           this->toStringWithPrecision(this->Variance()) + ")";
}

template <typename RealType>
void NormalRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Normal distribution: scale should be positive, but it's equal to " +
                                    std::to_string(scale));
    sigma = scale;
    StableDistribution<RealType>::SetScale(sigma * randlib::consts::inv_sqrt2);
}

template <typename RealType>
void NormalRand<RealType>::SetVariance(double variance)
{
    if (!(variance > 0.0))
        throw std::invalid_argument("Variance of Normal distribution should be positive, but it's equal to " +
                                    std::to_string(variance));
    SetScale(std::sqrt(variance));
}

template <typename RealType>
DoublePair NormalRand<RealType>::SufficientStatistic(RealType x) const
{
    return {x, x * x};
}

template <typename RealType>
DoublePair NormalRand<RealType>::SourceParameters() const
{
    return {this->mu, this->sigma * this->sigma};
}

template <typename RealType>
DoublePair NormalRand<RealType>::SourceToNatural(DoublePair sourceParameters) const
{
    double mean = sourceParameters.first;
    double scaleSq = sourceParameters.second;
    return {mean / scaleSq, -0.5 / scaleSq};
}

template <typename RealType>
double NormalRand<RealType>::LogNormalizer(DoublePair theta) const
{
    double theta1 = theta.first, theta2 = theta.second;
    double F = -0.25 * theta1 * theta1 / theta2;
    F += 0.5 * std::log(-randlib::consts::pi / theta2);
    return F;
}

template <typename RealType>
DoublePair NormalRand<RealType>::LogNormalizerGradient(DoublePair theta) const
{
    double theta1 = theta.first, theta2 = theta.second;
    double grad1 = -0.5 * theta1 / theta2;
    double grad2 = -0.5 / theta2 + grad1 * grad1;
    return {grad1, grad2};
}

template <typename RealType>
double NormalRand<RealType>::CarrierMeasure(RealType) const
{
    return 0.0;
}

template <typename RealType>
double NormalRand<RealType>::f(const RealType &x) const
{
    return this->pdfNormal(x);
}

template <typename RealType>
double NormalRand<RealType>::logf(const RealType &x) const
{
    return this->logpdfNormal(x);
}

template <typename RealType>
double NormalRand<RealType>::F(const RealType &x) const
{
    return this->cdfNormal(x);
}

template <typename RealType>
double NormalRand<RealType>::S(const RealType &x) const
{
    return this->cdfNormalCompl(x);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType NormalRand<RealType>::operator()(G &g) const
{
    return this->mu + sigma * StandardVariate(g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType NormalRand<RealType>::StandardVariate(G &g)
{
    return detail::standardNormal<RealType>(g);
}

template <typename RealType>
std::complex<double> NormalRand<RealType>::CFImpl(double t) const
{
    return this->cfNormal(t);
}

template <typename RealType>
RealType NormalRand<RealType>::quantileImpl(double p) const
{
    return this->quantileNormal(p);
}

template <typename RealType>
RealType NormalRand<RealType>::quantileImpl1m(double p) const
{
    return this->quantileNormal1m(p);
}

template <typename RealType>
long double NormalRand<RealType>::Moment(size_t n) const
{
    /// E[X^n] = Σ_{even k} C(n, k) μ^(n - k) σ^k (k - 1)!!
    long double sum = 0.0;
    for (size_t k = 0; k <= n; k += 2) {
        long double centralMoment = (k == 0) ? 1.0 : std::exp(k * std::log(sigma) + math::ldfact(k - 1));
        sum += math::binom(n, k) * std::pow(static_cast<long double>(this->mu), static_cast<long double>(n - k)) *
               centralMoment;
    }
    return sum;
}

template <typename RealType>
void NormalRand<RealType>::FitLocation(std::span<const RealType> sample)
{
    this->SetLocation(this->GetSampleMean(sample));
}

template <typename RealType>
void NormalRand<RealType>::FitLocation(std::span<const RealType> sample, DoublePair &confidenceInterval,
                                       double significanceLevel)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw this->fitError(FitErrorCode::WrongLevel,
                             "Input level is equal to " + this->toStringWithPrecision(significanceLevel));

    FitLocation(sample);

    int n = sample.size();
    NormalRand<RealType> NormalRV(0, 1);
    double halfAlpha = 0.5 * significanceLevel;
    double interval = NormalRV.Quantile1m(halfAlpha) * sigma / std::sqrt(n);
    confidenceInterval.first = this->mu - interval;
    confidenceInterval.second = this->mu + interval;
}

template <typename RealType>
void NormalRand<RealType>::FitVariance(std::span<const RealType> sample)
{
    this->SetVariance(this->GetSampleVariance(sample, this->mu));
}

template <typename RealType>
void NormalRand<RealType>::FitVariance(std::span<const RealType> sample, DoublePair &confidenceInterval,
                                       double significanceLevel, bool unbiased)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw this->fitError(FitErrorCode::WrongLevel,
                             "Input level is equal to " + this->toStringWithPrecision(significanceLevel));

    FitVariance(sample);

    size_t n = sample.size();
    double halfAlpha = 0.5 * significanceLevel;
    ChiSquaredRand ChiSqRV(n);
    double numerator = sigma * sigma * (unbiased ? n : (n - 1));
    confidenceInterval.first = numerator / ChiSqRV.Quantile1m(halfAlpha);
    confidenceInterval.second = numerator / ChiSqRV.Quantile(halfAlpha);
}

template <typename RealType>
void NormalRand<RealType>::FitScale(std::span<const RealType> sample, bool unbiased)
{
    if (unbiased == false)
        return FitVariance(sample);
    size_t n = sample.size();
    double halfN = 0.5 * n;
    double s = this->GetSampleVariance(sample, this->mu);
    s *= halfN;
    s = 0.5 * std::log(s);
    s -= std::lgammal(halfN + 0.5);
    s += std::lgammal(halfN);
    SetScale(std::exp(s));
}

template <typename RealType>
void NormalRand<RealType>::Fit(std::span<const RealType> sample, bool unbiased)
{
    double adjustment = 1.0;
    if (unbiased == true) {
        size_t n = sample.size();
        if (n <= 1)
            throw this->fitError(FitErrorCode::TooFewElements, "There should be at least 2 elements");
        adjustment = static_cast<double>(n) / (n - 1);
    }
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    this->SetLocation(stats.first);
    this->SetVariance(stats.second * adjustment);
}

template <typename RealType>
void NormalRand<RealType>::Fit(std::span<const RealType> sample, DoublePair &confidenceIntervalForMean,
                               DoublePair &confidenceIntervalForVariance, double significanceLevel, bool unbiased)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw this->fitError(FitErrorCode::WrongLevel,
                             "Input level is equal to " + this->toStringWithPrecision(significanceLevel));

    Fit(sample, unbiased);

    size_t n = sample.size();
    double sigmaAdj = unbiased ? sigma : (sigma * n) / (n - 1);
    /// calculate confidence interval for mean
    double halfAlpha = 0.5 * significanceLevel;
    StudentTRand<RealType> tRV(n - 1);
    double interval = tRV.Quantile1m(halfAlpha) * sigmaAdj / std::sqrt(n);
    confidenceIntervalForMean.first = this->mu - interval;
    confidenceIntervalForMean.second = this->mu + interval;

    /// calculate confidence interval for variance
    ChiSquaredRand<RealType> ChiSqRV(n - 1);
    double numerator = (n - 1) * sigmaAdj * sigmaAdj;
    confidenceIntervalForVariance.first = numerator / ChiSqRV.Quantile1m(halfAlpha);
    confidenceIntervalForVariance.second = numerator / ChiSqRV.Quantile(halfAlpha);
}

template <typename RealType>
NormalRand<RealType> NormalRand<RealType>::FitLocationBayes(std::span<const RealType> sample,
                                                            const NormalRand<RealType> &priorDistribution, bool MAP)
{
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = GetPrecision();
    double numerator = this->GetSampleSum(sample) * tau + tau0 * mu0;
    double denominator = sample.size() * tau + tau0;
    NormalRand<RealType> posteriorDistribution(numerator / denominator, 1.0 / denominator);
    this->SetLocation(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template <typename RealType>
InverseGammaRand<RealType> NormalRand<RealType>::FitVarianceBayes(std::span<const RealType> sample,
                                                                  const InverseGammaRand<RealType> &priorDistribution,
                                                                  bool MAP)
{
    double halfN = 0.5 * sample.size();
    double alphaPrior = priorDistribution.GetShape();
    double betaPrior = priorDistribution.GetRate();
    double alphaPosterior = alphaPrior + halfN;
    double betaPosterior = betaPrior + halfN * this->GetSampleVariance(sample, this->mu);
    InverseGammaRand<RealType> posteriorDistribution(alphaPosterior, betaPosterior);
    SetVariance(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template <typename RealType>
NormalInverseGammaRand<RealType>
NormalRand<RealType>::FitBayes(std::span<const RealType> sample,
                               const NormalInverseGammaRand<RealType> &priorDistribution, bool MAP)
{
    size_t n = sample.size();
    double alphaPrior = priorDistribution.GetShape();
    double betaPrior = priorDistribution.GetRate();
    double muPrior = priorDistribution.GetLocation();
    double lambdaPrior = priorDistribution.GetPrecision();
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    double lambdaPosterior = lambdaPrior + n;
    double muPosterior = (lambdaPrior * muPrior + n * stats.first) / lambdaPosterior;
    double halfN = 0.5 * n;
    double alphaPosterior = alphaPrior + halfN;
    double aux = muPrior - stats.first;
    double betaPosterior = betaPrior + halfN * (stats.second + lambdaPrior / lambdaPosterior * aux * aux);
    NormalInverseGammaRand<RealType> posteriorDistribution(muPosterior, lambdaPosterior, alphaPosterior, betaPosterior);
    DoublePair newParams = MAP ? static_cast<DoublePair>(posteriorDistribution.Mode())
                               : static_cast<DoublePair>(posteriorDistribution.Mean());
    this->SetLocation(newParams.first);
    this->SetVariance(newParams.second);
    return posteriorDistribution;
}

} // namespace randlib

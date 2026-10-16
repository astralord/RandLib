#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>
#include <randlib/core/exponential_family.hpp>

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN

namespace randlib
{

/**
 * @brief The PoissonRand class <BR>
 * Poisson distribution
 *
 * P(X = k) = λ^k * exp(-λ) / k!
 *
 * Notation: X ~ Po(λ)
 */
template <typename IntType>
class PoissonRand final : public DiscreteDistribution<IntType>, public ExponentialFamily<IntType, double>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<IntType> outputData) const override;

    double lambda = 1; ///< rate λ
    double logLambda = 0; ///< ln(λ)
    double Fmu = 2 * randlib::consts::inv_e; ///< P(X < [λ])
    double Pmu = randlib::consts::inv_e; ///< P(X = [λ])

    double mu = 1, delta = 6;
    double zeta{};
    long double c1{}, c2{}, c3{}, c4{}, c{};
    double sqrtMu = 1, sqrtMupHalfDelta = 2;
    double lfactMu = 0;

public:
    explicit PoissonRand(double rate = 1.0);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] IntType MinValue() const override { return 0; }
    [[nodiscard]] IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

private:
    void SetGeneratorConstants();

public:
    void SetRate(double rate);
    [[nodiscard]] inline double GetRate() const noexcept { return lambda; }

    double SufficientStatistic(IntType x) const override;
    double SourceParameters() const override;
    double SourceToNatural(double sourceParameters) const override;
    double NaturalParameters() const override;
    double LogNormalizer(double theta) const override;
    double LogNormalizerGradient(double theta) const override;
    double CarrierMeasure(IntType x) const override;
    double CrossEntropyAdjusted(double parameters) const override;
    double EntropyAdjusted() const override;

    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;
    [[nodiscard]] double S(const IntType &k) const override;

private:
    double acceptanceFunction(IntType X) const;
    bool generateByInversion() const;
    template <std::uniform_random_bit_generator G>
    IntType variateRejection(G &g) const;
    template <std::uniform_random_bit_generator G>
    IntType variateInversion(G &g) const;

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static IntType Variate(double rate, G &g);

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Median() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @fn Fit
     * fit rate λ via maximum-likelihood method
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void Fit(std::span<const IntType> sample);
    /**
     * @brief Fit
     * @param sample
     * @param confidenceInterval
     * @param significanceLevel
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void Fit(std::span<const IntType> sample, DoublePair &confidenceInterval,
                                           double significanceLevel);
    /**
     * @fn FitBayes
     * fit rate λ via Bayes estimation
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior Gamma distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER GammaRand<> FitBayes(std::span<const IntType> sample,
                                                       const GammaDistribution<> &priorDistribution, bool MAP = false);
};

} // namespace randlib

#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/exponential.hpp>

namespace randlib
{

template <typename IntType>
PoissonRand<IntType>::PoissonRand(double rate)
{
    SetRate(rate);
}

template <typename IntType>
std::string PoissonRand<IntType>::Name() const
{
    return "Poisson(" + this->toStringWithPrecision(GetRate()) + ")";
}

template <typename IntType>
void PoissonRand<IntType>::SetGeneratorConstants()
{
    delta = std::round(std::sqrt(2 * mu * std::log(randlib::consts::inv_pi * 128 * mu)));
    delta = std::max(6.0, std::min(mu, delta));

    c1 = std::sqrt(0.5 * randlib::consts::pi * mu * randlib::consts::e);

    c2 = 1.0 / (2 * mu + delta);
    c2 += 0.5 * (randlib::consts::lnpi - randlib::consts::ln2 + std::log(mu + 0.5 * delta));
    c2 = c1 + std::exp(c2);

    c3 = c2 + 1.0;
    c4 = c3 + 1.0;

    zeta = (4.0 * mu) / delta + 2;
    c = std::exp(-(2.0 + delta) / zeta);
    c *= zeta;
    c += c4;

    sqrtMu = std::sqrt(mu);
    sqrtMupHalfDelta = std::sqrt(mu + 0.5 * delta);
    lfactMu = math::lfact(mu);
}

template <typename IntType>
void PoissonRand<IntType>::SetRate(double rate)
{
    if (!(rate > 0.0))
        throw std::invalid_argument("Poisson distribution: rate should be positive, but it's equal to " +
                                    std::to_string(rate));
    lambda = rate;

    logLambda = std::log(lambda);
    mu = std::floor(lambda);
    Fmu = F(mu);
    Pmu = this->P(mu);

    if (!generateByInversion())
        SetGeneratorConstants();
}

template <typename IntType>
double PoissonRand<IntType>::SufficientStatistic(IntType x) const
{
    return x;
}

template <typename IntType>
double PoissonRand<IntType>::SourceParameters() const
{
    return lambda;
}

template <typename IntType>
double PoissonRand<IntType>::SourceToNatural(double sourceParameters) const
{
    return std::log(sourceParameters);
}

template <typename IntType>
double PoissonRand<IntType>::NaturalParameters() const
{
    return logLambda;
}

template <typename IntType>
double PoissonRand<IntType>::LogNormalizer(double theta) const
{
    return std::exp(theta);
}

template <typename IntType>
double PoissonRand<IntType>::LogNormalizerGradient(double theta) const
{
    return std::exp(theta);
}

template <typename IntType>
double PoissonRand<IntType>::CarrierMeasure(IntType x) const
{
    return -math::lfact(x);
}

template <typename IntType>
double PoissonRand<IntType>::CrossEntropyAdjusted(double parameters) const
{
    return parameters - lambda * std::log(parameters);
}

template <typename IntType>
double PoissonRand<IntType>::EntropyAdjusted() const
{
    return lambda - lambda * logLambda;
}

template <typename IntType>
double PoissonRand<IntType>::logP(const IntType &k) const
{
    if (k < 0)
        return -INFINITY;
    if (k == 0)
        return -lambda;
    /// saddle-point form; k log λ - λ - log(k!) cancels catastrophically for large k
    return -math::stirlingError(k) - math::bd0(k, lambda) - 0.5 * std::log(2 * randlib::consts::pi * k);
}

template <typename IntType>
double PoissonRand<IntType>::F(const IntType &k) const
{
    return (k >= 0.0) ? math::qgamma(k + 1, lambda, logLambda) : 0.0;
}

template <typename IntType>
double PoissonRand<IntType>::S(const IntType &k) const
{
    return (k >= 0.0) ? math::pgamma(k + 1, lambda, logLambda) : 1.0;
}

template <typename IntType>
double PoissonRand<IntType>::acceptanceFunction(IntType X) const
{
    if (X == 0)
        return 0.0;
    double q = X * logLambda;
    q += lfactMu;
    q -= math::lfact(X + mu);
    return q;
}

template <typename IntType>
bool PoissonRand<IntType>::generateByInversion() const
{
    /// the inversion generator is much faster than rejection,
    /// however precision loss for large rate increases drastically
    return lambda < 10;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType PoissonRand<IntType>::variateRejection(G &g) const
{
    size_t iter = 0;
    IntType X = 0;
    do {
        bool reject = false;
        float W = 0.0;
        float U = c * detail::uniform01<float>(g);
        if (U <= c1) {
            float N = detail::standardNormal<float>(g);
            float Y = -std::fabs(N) * sqrtMu;
            X = std::floor(Y);
            if (X < -mu) {
                reject = true;
            } else {
                W = -0.5 * (N * N - 1.0);
            }
        } else if (U <= c2) {
            float N = detail::standardNormal<float>(g);
            float Y = 1.0 + std::fabs(N) * sqrtMupHalfDelta;
            X = std::ceil(Y);
            if (X > delta) {
                reject = true;
            } else {
                W = Y * (2.0 - Y) / (2.0 * mu + delta);
            }
        } else if (U <= c3) {
            return mu;
        } else if (U <= c4) {
            X = 1;
        } else {
            float V = detail::standardExponential<float>(g);
            float Y = delta + V * zeta;
            X = std::ceil(Y);
            W = -(2.0 + Y) / zeta;
        }

        if (!reject && W - detail::standardExponential<float>(g) <= acceptanceFunction(X)) {
            return X + mu;
        }

    } while (++iter < ProbabilityDistribution<IntType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Poisson distribution: sampling failed");
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType PoissonRand<IntType>::variateInversion(G &g) const
{
    double U = detail::uniform01<double>(g);
    IntType k = mu;
    double s = Fmu, p = Pmu;
    if (s < U) {
        do {
            ++k;
            p *= lambda / k;
            s += p;
        } while (s < U && p > 0);
    } else {
        s -= p;
        while (k > 0 && s > U) {
            p /= lambda / k;
            --k;
            s -= p;
        }
    }
    return k;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType PoissonRand<IntType>::operator()(G &g) const
{
    return generateByInversion() ? variateInversion(g) : variateRejection(g);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType PoissonRand<IntType>::Variate(double rate, G &g)
{
    /// check validness of parameter
    if (!(rate > 0.0))
        throw std::invalid_argument("Poisson distribution: rate should be positive");
    if (rate > 1000) {
        /// approximate with normal distribution
        float X = detail::standardNormal<float>(g);
        return std::floor(rate + std::sqrt(rate) * X);
    }
    int k = -1;
    double s = 0;
    do {
        s += detail::standardExponential<double>(g);
        ++k;
    } while (s < rate);
    return k;
}

template <typename IntType>
void PoissonRand<IntType>::sample(detail::EngineRef g, std::span<IntType> outputData) const
{
    if (generateByInversion()) {
        for (IntType &var : outputData)
            var = variateInversion(g);
    } else {
        for (IntType &var : outputData)
            var = variateRejection(g);
    }
}

template <typename IntType>
long double PoissonRand<IntType>::Mean() const
{
    return lambda;
}

template <typename IntType>
long double PoissonRand<IntType>::Variance() const
{
    return lambda;
}

template <typename IntType>
std::complex<double> PoissonRand<IntType>::CFImpl(double t) const
{
    std::complex<double> y(std::cos(t) - 1.0, std::sin(t));
    return std::exp(lambda * y);
}

template <typename IntType>
IntType PoissonRand<IntType>::Median() const
{
    /// this value is approximate
    return std::max(std::floor(lambda + 1.0 / 3 - 0.02 / lambda), 0.0);
}

template <typename IntType>
IntType PoissonRand<IntType>::Mode() const
{
    if (math::detail::areClose(mu, lambda)) {
        return (Pmu < this->P(mu + 1)) ? mu + 1 : mu;
    }
    return mu;
}

template <typename IntType>
long double PoissonRand<IntType>::Skewness() const
{
    return 1.0 / std::sqrt(lambda);
}

template <typename IntType>
long double PoissonRand<IntType>::ExcessKurtosis() const
{
    return 1.0 / lambda;
}

template <typename IntType>
void PoissonRand<IntType>::Fit(std::span<const IntType> sample)
{
    if (!this->allElementsAreNonNegative(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->NON_NEGATIVITY_VIOLATION);
    SetRate(this->GetSampleMean(sample));
}

template <typename IntType>
void PoissonRand<IntType>::Fit(std::span<const IntType> sample, DoublePair &confidenceInterval,
                               double significanceLevel)
{
    size_t n = sample.size();

    if (significanceLevel <= 0 || significanceLevel > 1)
        throw this->fitError(FitErrorCode::WrongLevel,
                             "Alpha is equal to " + this->toStringWithPrecision(significanceLevel));

    Fit(sample);

    double halfAlpha = 0.5 * significanceLevel;
    ErlangRand<double> ErlangRV(n);
    confidenceInterval.first = ErlangRV.Quantile(halfAlpha);
    ErlangRV.SetShape(n + 1);
    confidenceInterval.second = ErlangRV.Quantile1m(halfAlpha);
}

template <typename IntType>
GammaRand<> PoissonRand<IntType>::FitBayes(std::span<const IntType> sample,
                                           const GammaDistribution<> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNonNegative(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->NON_NEGATIVITY_VIOLATION);
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    GammaRand<> posteriorDistribution(alpha + this->GetSampleSum(sample), beta + sample.size());
    SetRate(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

} // namespace randlib

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END

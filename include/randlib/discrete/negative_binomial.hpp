#pragma once

#include <algorithm>
#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>
#include <randlib/continuous/gamma.hpp>

namespace randlib
{

/**
 * @brief The NegativeBinomialDistribution class <BR>
 * Abstract class for Negative binomial distribution
 *
 * P(X = k) = C(k + r - 1, k) p^r (1-p)^k
 *
 * Notation: X ~ NB(r, p)
 *
 * Related distributions: <BR>
 * If X ~ NB(1, p), then X ~ Geometric(p)
 * If Y ~ Γ(r, p / (1 - p), then Po(Y) ~ NB(r, p)
 */
template <typename IntType, typename T>
class NegativeBinomialDistribution : public DiscreteDistribution<IntType>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<IntType> outputData) const override;

protected:
    double p = 0.5; ///< probability of failure
    double q = 0.5; ///< probability of success
    double logProb = -randlib::consts::ln2; ///< log(p)
    double log1mProb = -randlib::consts::ln2; ///< log(q)

private:
    T r = 1; ///< number of failures until the experiment is stopped
    double pdfCoef = -randlib::consts::ln2; ///< coefficient for faster pdf calculation
    double qDivP = 1; ///< q / p
    static constexpr int tableSize = 16;
    double table[tableSize];
    GammaRand<float> GammaRV{};

protected:
    NegativeBinomialDistribution(T number = 1, double probability = 0.5);

public:
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] IntType MinValue() const override { return 0; }
    [[nodiscard]] IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

protected:
    void SetParameters(T number, double probability);

public:
    [[nodiscard]] inline double GetProbability() const noexcept { return p; }
    [[nodiscard]] inline T GetNumber() const noexcept { return r; }

    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;
    [[nodiscard]] double S(const IntType &k) const override;

protected:
    enum GENERATOR_ID
    {
        TABLE,
        EXPONENTIAL,
        GAMMA_POISSON
    };

    /**
     * If r is small, we use two different generators for two different cases:
     * If p < 0.08 then the tail is too heavy (probability to be in main body is smaller than 0.75),
     * then we return highest integer, smaller than variate from exponential distribution.
     * Otherwise we choose table method
     * @return id of generator
     */
    [[nodiscard]] GENERATOR_ID GetIdOfUsedGenerator() const
    {
        if ((r < 10 || GammaRV.Mean() > 10) && r == std::round(r))
            return (p < 0.08) ? EXPONENTIAL : TABLE;
        return GAMMA_POISSON;
    }

    template <std::uniform_random_bit_generator G>
    IntType variateGeometricByTable(G &g) const;
    template <std::uniform_random_bit_generator G>
    IntType variateGeometricThroughExponential(G &g) const;

private:
    template <std::uniform_random_bit_generator G>
    IntType variateByTable(G &g) const;
    template <std::uniform_random_bit_generator G>
    IntType variateThroughExponential(G &g) const;
    template <std::uniform_random_bit_generator G>
    IntType variateThroughGammaPoisson(G &g) const;

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER BetaRand<>
    FitProbabilityBayes(std::span<const IntType> sample, const BetaDistribution<> &priorDistribution, bool MAP = false);
};

/**
 * @brief The NegativeBinomialRand class <BR>
 * Negative binomial distribution
 */
template <typename IntType, typename T>
class NegativeBinomialRand : public NegativeBinomialDistribution<IntType, T>
{
public:
    NegativeBinomialRand(T number = 1, double probability = 0.5)
        : NegativeBinomialDistribution<IntType, T>(number, probability)
    {
    }
    [[nodiscard]] std::string Name() const override;

    using NegativeBinomialDistribution<IntType, T>::SetParameters;

    static constexpr char TOO_SMALL_VARIANCE[] = "Sample variance should be greater than sample mean";

    /**
     * set number and probability, estimated via maximum-likelihood method
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void Fit(std::span<const IntType> sample);
};

template <typename IntType = int>
using PascalRand = NegativeBinomialRand<IntType, int>;

template <typename IntType = int>
using PolyaRand = NegativeBinomialRand<IntType, double>;

} // namespace randlib

#include <randlib/continuous/beta.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/discrete/poisson.hpp>

namespace randlib
{

template <typename IntType, typename T>
NegativeBinomialDistribution<IntType, T>::NegativeBinomialDistribution(T number, double probability)
{
    SetParameters(number, probability);
}

template <typename IntType, typename T>
void NegativeBinomialDistribution<IntType, T>::SetParameters(T number, double probability)
{
    if (!(number > 0))
        throw std::invalid_argument("Negative-Binomial distribution: number parameter should be positive");
    if (!(probability > 0.0 && probability <= 1.0))
        throw std::invalid_argument(
            "Negative-Binomial distribution: probability parameter should be in interval (0, 1]");
    r = number;
    p = probability;
    q = 1.0 - p;
    logProb = std::log(p);
    log1mProb = std::log1pl(-p);
    if (p == 1.0) {
        qDivP = 0.0;
        pdfCoef = -std::lgamma(r);
        return;
    }
    GammaRV.SetParameters(r, p / q);
    qDivP = GammaRV.GetScale();
    pdfCoef = r * logProb;
    pdfCoef -= GammaRV.GetLogGammaShape();

    if (GetIdOfUsedGenerator() == TABLE) {
        /// table method
        table[0] = p;
        double prod = p;
        for (int i = 1; i < tableSize; ++i) {
            prod *= q;
            table[i] = table[i - 1] + prod;
        }
    }
}

template <typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::P(const IntType &k) const
{
    return (k < 0) ? 0.0 : std::exp(logP(k));
}

template <typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::logP(const IntType &k) const
{
    if (k < 0)
        return -INFINITY;
    if (p == 1.0)
        return (k == 0) ? 0.0 : -INFINITY;
    double y = std::lgammal(r + k);
    y -= math::lfact(k);
    y += k * log1mProb;
    y += pdfCoef;
    return y;
}

template <typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::F(const IntType &k) const
{
    if (k < 0)
        return 0.0;
    if (p == 1.0)
        return 1.0;
    int kp1 = k + 1;
    double logBetaFun = math::logBeta(r, kp1);
    return math::ibeta(p, r, kp1, logBetaFun, logProb, log1mProb);
}

template <typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::S(const IntType &k) const
{
    if (k < 0)
        return 1.0;
    if (p == 1.0)
        return 0.0;
    int kp1 = k + 1;
    double logBetaFun = math::logBeta(kp1, r);
    return math::ibeta(q, kp1, r, logBetaFun, log1mProb, logProb);
}

template <typename IntType, typename T>
template <std::uniform_random_bit_generator G>
IntType NegativeBinomialDistribution<IntType, T>::variateThroughGammaPoisson(G &g) const
{
    /// For small shapes the gamma rate underflows to zero, the limit in which the Poisson variate is zero
    const float rate = GammaRV(g);
    return (rate > 0) ? PoissonRand<IntType>::Variate(rate, g) : 0;
}

template <typename IntType, typename T>
template <std::uniform_random_bit_generator G>
IntType NegativeBinomialDistribution<IntType, T>::variateGeometricByTable(G &g) const
{
    double U = detail::uniform01<double>(g);
    /// handle tail by recursion
    if (U > table[tableSize - 1])
        return tableSize + variateGeometricByTable(g);
    /// handle the main body
    IntType x = 0;
    while (U > table[x])
        ++x;
    return x;
}

template <typename IntType, typename T>
template <std::uniform_random_bit_generator G>
IntType NegativeBinomialDistribution<IntType, T>::variateGeometricThroughExponential(G &g) const
{
    float X = -detail::standardExponential<float>(g) / log1mProb;
    if (!std::isfinite(X) || X >= static_cast<long double>(std::numeric_limits<IntType>::max()))
        return std::numeric_limits<IntType>::max();
    return static_cast<IntType>(X);
}

template <typename IntType, typename T>
template <std::uniform_random_bit_generator G>
IntType NegativeBinomialDistribution<IntType, T>::variateByTable(G &g) const
{
    IntType var = 0;
    for (int i = 0; i < r; ++i) {
        IntType add = variateGeometricByTable(g);
        if (add > std::numeric_limits<IntType>::max() - var)
            return std::numeric_limits<IntType>::max();
        var += add;
    }
    return var;
}

template <typename IntType, typename T>
template <std::uniform_random_bit_generator G>
IntType NegativeBinomialDistribution<IntType, T>::variateThroughExponential(G &g) const
{
    IntType var = 0;
    for (int i = 0; i < r; ++i) {
        IntType add = variateGeometricThroughExponential(g);
        if (add > std::numeric_limits<IntType>::max() - var)
            return std::numeric_limits<IntType>::max();
        var += add;
    }
    return var;
}

template <typename IntType, typename T>
template <std::uniform_random_bit_generator G>
IntType NegativeBinomialDistribution<IntType, T>::operator()(G &g) const
{
    if (p == 1.0)
        return 0;
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == TABLE)
        return variateByTable(g);
    return (genId == EXPONENTIAL) ? variateThroughExponential(g) : variateThroughGammaPoisson(g);
}

template <typename IntType, typename T>
void NegativeBinomialDistribution<IntType, T>::sample(detail::EngineRef g, std::span<IntType> outputData) const
{
    if (p == 1.0) {
        std::fill(outputData.begin(), outputData.end(), 0);
        return;
    }
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == TABLE) {
        for (IntType &var : outputData)
            var = variateByTable(g);
    } else if (genId == EXPONENTIAL) {
        for (IntType &var : outputData)
            var = variateThroughExponential(g);
    } else {
        for (IntType &var : outputData)
            var = variateThroughGammaPoisson(g);
    }
}

template <typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::Mean() const
{
    return qDivP * r;
}

template <typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::Variance() const
{
    return qDivP * r / p;
}

template <typename IntType, typename T>
std::complex<double> NegativeBinomialDistribution<IntType, T>::CFImpl(double t) const
{
    std::complex<double> denominator(1.0 - q * std::cos(t), -q * std::sin(t));
    return std::pow(p / denominator, r);
}

template <typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::Mode() const
{
    return (r > 1) ? std::floor((r - 1) * qDivP) : 0;
}

template <typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::Skewness() const
{
    if (p == 1.0)
        return NAN;
    return (1 + q) / std::sqrt(q * r);
}

template <typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::ExcessKurtosis() const
{
    if (p == 1.0)
        return NAN;
    long double kurtosis = p / qDivP;
    kurtosis += 6;
    return kurtosis / r;
}

template <typename IntType, typename T>
BetaRand<> NegativeBinomialDistribution<IntType, T>::FitProbabilityBayes(std::span<const IntType> sample,
                                                                         const BetaDistribution<> &priorDistribution,
                                                                         bool MAP)
{
    if (!this->allElementsAreNonNegative(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->NON_NEGATIVITY_VIOLATION);
    int n = sample.size();
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    BetaRand<> posteriorDistribution(alpha + r * n, beta + this->GetSampleSum(sample));
    SetParameters(r, MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}
template <typename IntType, typename T>
std::string NegativeBinomialRand<IntType, T>::Name() const
{
    if (std::is_integral_v<T>)
        return "Pascal(" + this->toStringWithPrecision(this->GetNumber()) + ", " +
               this->toStringWithPrecision(this->GetProbability()) + ")";
    return "Polya(" + this->toStringWithPrecision(this->GetNumber()) + ", " +
           this->toStringWithPrecision(this->GetProbability()) + ")";
}

template <typename IntType, typename T>
constexpr char NegativeBinomialRand<IntType, T>::TOO_SMALL_VARIANCE[];

template <typename IntType, typename T>
void NegativeBinomialRand<IntType, T>::Fit(std::span<const IntType> sample)
{
    /// Check positivity of sample
    if (!this->allElementsAreNonNegative(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->NON_NEGATIVITY_VIOLATION);
    /// Initial guess by method of moments
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    double mean = stats.first, variance = stats.second;
    /// Method can't be applied in case of variance smaller than mean
    if (variance <= mean)
        throw this->fitError(FitErrorCode::NotApplicable, this->TOO_SMALL_VARIANCE);
    double guess = mean * mean / (variance - mean);
    size_t n = sample.size();
    auto score = [sample, mean, n](double x) {
        double first = 0.0, second = 0.0;
        for (const IntType &var : sample) {
            first += math::digamma(var + x);
            second += math::trigamma(var + x);
        }
        first -= n * (math::digammamLog(x) + std::log(x + mean));
        second -= n * (math::trigamma(x) - mean / (x * (mean + x)));
        return DoublePair(first, second);
    };
    double logGuess = std::log(guess);
    if (!math::findRootNewtonFirstOrder<double>(
            [&score](double logNumber) {
                double number = std::exp(logNumber);
                DoublePair value = score(number);
                value.second *= number;
                return value;
            },
            logGuess))
        throw this->fitError(FitErrorCode::Undefined, "Error in root-finding algorithm");
    guess = std::exp(logGuess);
    if (!std::isfinite(guess) || guess <= 0.0)
        throw this->fitError(FitErrorCode::WrongReturn,
                             "Number should be positive, but returned value is " + this->toStringWithPrecision(guess));
    SetParameters(guess, guess / (guess + mean));
}

} // namespace randlib

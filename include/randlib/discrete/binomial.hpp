#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>
#include <randlib/core/exponential_family.hpp>
#include <randlib/discrete/geometric.hpp>

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_BEGIN

namespace randlib
{

/**
 * @brief The BinomialDistribution class <BR>
 * Abstract class for Binomial distribution
 *
 * Notation: X ~ Bin(n, p)
 *
 * Related distributions: <BR>
 * If X ~ Bin(1, p), then X ~ Bernoulli(p) <BR>
 * X ~ Multin(n, 1 - p, p)
 */
template <typename IntType>
class BinomialDistribution : public DiscreteDistribution<IntType>, public ExponentialFamily<IntType, double>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<IntType> outputData) const override;

protected:
    double p = 0.5; ///< probability of success
    double q = 0.5; ///< probability of failure
    double logProb = -randlib::consts::ln2; ///< log(p)
    double log1mProb = -randlib::consts::ln2; ///< log(q)

private:
    IntType n = 1; ///< number of experiments
    double np = 0.5; ///< n * p
    double lfactn = 0; ///< log(n!)

    double delta1{}, delta2{};
    double sigma1{}, sigma2{}, c{};
    double a1{}, a2{}, a3{}, a4{};
    double coefa3{}, coefa4{};

    double minpq = 0.5; ///< min(p, q)
    double pFloor = 0; ///< [n * min(p, q)] / n
    double logPFloor = -INFINITY; ///< log(pFloor)
    double logQFloor = 0; ///< log(1 - pFloor)
    double pRes = 0.5; ///< min(p, q) - pFloor
    double npFloor = 0; ///< [n * min(p, q)]
    double nqFloor = 0; ///< [n * max(p, q)]
    double logPnpInv = 0; ///< log(P([npFloor)) if p = pFloor

    GeometricRand<IntType> geometric{};

protected:
    BinomialDistribution(IntType number, double probability);

public:
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] IntType MinValue() const override { return 0; }
    [[nodiscard]] IntType MaxValue() const override { return n; }

private:
    void SetGeneratorConstants();

protected:
    void SetParameters(IntType number, double probability);

public:
    [[nodiscard]] inline IntType GetNumber() const noexcept { return n; }
    [[nodiscard]] inline double GetProbability() const noexcept { return p; }

    double SufficientStatistic(IntType x) const override;
    double SourceParameters() const override;
    double SourceToNatural(double sourceParameters) const override;
    double NaturalParameters() const override;
    double LogNormalizer(double theta) const override;
    double LogNormalizerGradient(double theta) const override;
    double CarrierMeasure(IntType x) const override;
    double EntropyAdjusted() const override;

private:
    /**
     * @param k
     * @return logarithm of probability to get k if p = pFloor
     */
    double logProbFloor(int k) const;

public:
    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;
    [[nodiscard]] double S(const IntType &k) const override;

private:
    enum GENERATOR_ID
    {
        BERNOULLI_SUM,
        WAITING,
        REJECTION
    };

    [[nodiscard]] GENERATOR_ID GetIdOfUsedGenerator() const
    {
        /// if (n is tiny and minpq is big) or (p ~= 0.5 and n is not that large),
        /// we just sum Bernoulli random variables
        if ((n <= 3) || (n <= 13 && minpq > 0.025 * (n + 6)) || (n <= 200 && math::detail::areClose(p, 0.5)))
            return BERNOULLI_SUM;

        /// for small [np] we use simple waiting algorithm
        if ((npFloor <= 12) || (pRes > 0 && npFloor <= 16))
            return WAITING;

        /// otherwise
        return REJECTION;
    }

    template <std::uniform_random_bit_generator G>
    IntType variateRejection(G &g) const;
    template <std::uniform_random_bit_generator G>
    IntType variateWaiting(IntType number, G &g) const;
    template <std::uniform_random_bit_generator G>
    static IntType variateWaiting(IntType number, double probability, G &g);
    template <std::uniform_random_bit_generator G>
    static IntType variateBernoulliSum(IntType number, double probability, G &g);

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static IntType Variate(IntType number, double probability, G &g);

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Median() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

    /**
     * @return log(n!)
     */
    [[nodiscard]] inline double GetLogFactorialN() const noexcept { return lfactn; }
    /**
     * @return log(p)
     */
    [[nodiscard]] inline double GetLogProbability() const noexcept { return logProb; }
    /**
     * @return log(1-p)
     */
    [[nodiscard]] inline double GetLog1mProbability() const noexcept { return log1mProb; }

private:
    std::complex<double> CFImpl(double t) const override;

public:
    /**
     * Fit probability p with maximum-likelihood estimation
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void FitProbability(std::span<const IntType> sample);

    /**
     * Fit probability p with prior assumption p ~ Beta(α, β)
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER BetaRand<>
    FitProbabilityBayes(std::span<const IntType> sample, const BetaDistribution<> &priorDistribution, bool MAP = false);

    /**
     * Fit probability p with minimax estimator
     * @param sample
     * @return posterior distribution
     */
    RANDLIB_DEPRECATED_FIT_MEMBER BetaRand<> FitProbabilityMinimax(std::span<const IntType> sample);
};

/**
 * @brief The BinomialRand class <BR>
 * Binomial distribution
 */
template <typename IntType>
class BinomialRand final : public BinomialDistribution<IntType>
{
public:
    BinomialRand(int number = 1, double probability = 0.5) : BinomialDistribution<IntType>(number, probability) {}
    [[nodiscard]] std::string Name() const override;
    using BinomialDistribution<IntType>::SetParameters;
};

} // namespace randlib

#include <randlib/continuous/beta.hpp>
#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/exponential.hpp>
#include <randlib/discrete/bernoulli.hpp>

namespace randlib
{

template <typename IntType>
BinomialDistribution<IntType>::BinomialDistribution(IntType number, double probability)
{
    SetParameters(number, probability);
}

template <typename IntType>
void BinomialDistribution<IntType>::SetGeneratorConstants()
{
    minpq = std::min(p, q);
    npFloor = std::floor(n * minpq);
    pFloor = npFloor / n;
    pRes = math::detail::areClose(npFloor, n * minpq) ? 0.0 : minpq - pFloor;

    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == BERNOULLI_SUM)
        return;
    else if (genId == WAITING) {
        geometric.SetProbability(minpq);
        return;
    }

    nqFloor = n - npFloor;
    double qFloor = 1.0 - pFloor;
    if (pRes > 0)
        geometric.SetProbability(pRes / qFloor);

    /// Set deltas
    double npq = npFloor * qFloor;
    double coef = 128.0 * n / randlib::consts::pi;
    delta1 = coef * pFloor / (81.0 * qFloor);
    delta1 = npq * std::log(delta1);
    if (delta1 > 1.0)
        delta1 = std::sqrt(delta1);
    else
        delta1 = 1.0;
    delta2 = coef * qFloor / pFloor;
    delta2 = npq * std::log(delta2);
    if (delta2 > 1.0)
        delta2 = std::sqrt(delta2);
    else
        delta2 = 1.0;

    /// Set sigmas and c
    double npqSqrt = std::sqrt(npq);
    sigma1 = npqSqrt * (1.0 + 0.25 * delta1 / npFloor);
    sigma2 = npqSqrt * (1.0 + 0.25 * delta2 / nqFloor);
    c = 2.0 * delta1 / npFloor;

    /// Set a's
    a1 = 0.5 * std::exp(c) * sigma1 * randlib::consts::sqrt2pi;

    a2 = 0.5 * sigma2 * randlib::consts::sqrt2pi;
    a2 += a1;

    coefa3 = 0.5 * delta1 / (sigma1 * sigma1);
    a3 = 1.0 / nqFloor - coefa3;
    a3 *= delta1;
    a3 = std::exp(a3);
    a3 /= coefa3;
    a3 += a2;

    coefa4 = 0.5 * delta2 / (sigma2 * sigma2);
    a4 = std::exp(-delta2 * coefa4) / coefa4;
    a4 += a3;

    logPFloor = std::log(pFloor);
    logQFloor = (pFloor == qFloor) ? logPFloor : std::log1pl(-pFloor);

    logPnpInv = logProbFloor(npFloor);
}

template <typename IntType>
void BinomialDistribution<IntType>::SetParameters(IntType number, double probability)
{
    if (!(probability >= 0.0 && probability <= 1.0))
        throw std::invalid_argument("Binomial distribution: probability parameter should in interval [0, 1]");
    if (!(number > 0))
        throw std::invalid_argument("Binomial distribution: number should be positive");
    n = number;
    p = probability;
    q = 1.0 - p;
    np = n * p;
    lfactn = math::lfact(n);
    logProb = std::log(p);
    log1mProb = std::log1pl(-p);
    SetGeneratorConstants();
}

template <typename IntType>
double BinomialDistribution<IntType>::SufficientStatistic(IntType x) const
{
    return x;
}

template <typename IntType>
double BinomialDistribution<IntType>::SourceParameters() const
{
    return this->p;
}

template <typename IntType>
double BinomialDistribution<IntType>::SourceToNatural(double sourceParameters) const
{
    return std::log(sourceParameters) - std::log1p(-sourceParameters);
}

template <typename IntType>
double BinomialDistribution<IntType>::NaturalParameters() const
{
    return this->logProb - this->log1mProb;
}

template <typename IntType>
double BinomialDistribution<IntType>::LogNormalizer(double theta) const
{
    double F = n * math::softplus(theta);
    F -= this->lfactn;
    return F;
}

template <typename IntType>
double BinomialDistribution<IntType>::LogNormalizerGradient(double theta) const
{
    double expTheta = std::exp(theta);
    return this->n * expTheta / (1.0 + expTheta);
}

template <typename IntType>
double BinomialDistribution<IntType>::CarrierMeasure(IntType x) const
{
    double k = -math::lfact(x);
    k -= math::lfact(n - x);
    return k;
}

template <typename IntType>
double BinomialDistribution<IntType>::EntropyAdjusted() const
{
    double H = this->p * this->logProb;
    H += this->q * this->log1mProb;
    H *= this->n;
    H += this->lfactn;
    return -H;
}

template <typename IntType>
double BinomialDistribution<IntType>::logProbFloor(int k) const
{
    double y = lfactn;
    y -= math::lfact(n - k);
    y -= math::lfact(k);
    y += k * logPFloor;
    y += (n - k) * logQFloor;
    return y;
}

template <typename IntType>
double BinomialDistribution<IntType>::P(const IntType &k) const
{
    return (k < 0 || k > n) ? 0.0 : std::exp(logP(k));
}

template <typename IntType>
double BinomialDistribution<IntType>::logP(const IntType &k) const
{
    if (k < 0 || k > n)
        return -INFINITY;
    if (k == 0)
        return n * log1mProb;
    if (k == n)
        return n * logProb;
    /// saddle-point form; the direct sum of logarithms cancels catastrophically for large n
    double y = math::stirlingError(n) - math::stirlingError(k) - math::stirlingError(n - k);
    y -= math::bd0(k, n * p) + math::bd0(n - k, n * q);
    return y - 0.5 * std::log(2 * randlib::consts::pi * k * (1.0 - static_cast<double>(k) / n));
}

template <typename IntType>
double BinomialDistribution<IntType>::F(const IntType &k) const
{
    if (k < 0)
        return 0.0;
    if (k >= n)
        return 1.0;
    int nmk = n - k, kp1 = k + 1;
    double logBetaFun = math::lfact(n - kp1);
    logBetaFun += math::lfact(k);
    logBetaFun -= lfactn;
    return math::ibeta(q, nmk, kp1, logBetaFun, log1mProb, logProb);
}

template <typename IntType>
double BinomialDistribution<IntType>::S(const IntType &k) const
{
    if (k < 0)
        return 1.0;
    if (k >= n)
        return 0.0;
    int nmk = n - k, kp1 = k + 1;
    double logBetaFun = math::logBeta(kp1, nmk);
    return math::ibeta(p, kp1, nmk, logBetaFun, logProb, log1mProb);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BinomialDistribution<IntType>::variateRejection(G &g) const
{
    /// a rejection algorithm by Devroye and Naderlsamanl (1980)
    /// p.533. Non-Uniform Random Variate Generation. Luc Devroye
    /// it can be used only when n * p is integer and p < 0.5
    bool reject = true;
    size_t iter = 0;
    float Y, V = 0;
    IntType X = 0;
    do {
        float U = a4 * detail::uniform01<float>(g);
        if (U <= a1) {
            float N = detail::standardNormal<float>(g);
            Y = sigma1 * std::fabs(N);
            reject = (Y >= delta1);
            if (!reject) {
                float W = detail::standardExponential<float>(g);
                X = std::floor(Y);
                V = -W - 0.5 * N * N + c;
            }
        } else if (U <= a2) {
            float N = detail::standardNormal<float>(g);
            Y = sigma2 * std::fabs(N);
            reject = (Y >= delta2);
            if (!reject) {
                float W = detail::standardExponential<float>(g);
                X = std::floor(-Y);
                V = -W - 0.5 * N * N;
            }
        } else if (U <= a3) {
            float W1 = detail::standardExponential<float>(g);
            float W2 = detail::standardExponential<float>(g);
            Y = delta1 + W1 / coefa3;
            X = std::floor(Y);
            V = -W2 - coefa3 * Y + delta1 / nqFloor;
            reject = false;
        } else {
            float W1 = detail::standardExponential<float>(g);
            float W2 = detail::standardExponential<float>(g);
            Y = delta2 + W1 / coefa4;
            X = std::floor(-Y);
            V = -W2 - coefa4 * Y;
            reject = false;
        }

        if (!reject) {
            X += npFloor;
            if (X >= 0 && X <= n && V <= logProbFloor(X) - logPnpInv)
                return X;
        }
    } while (++iter <= ProbabilityDistribution<IntType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Binomial distribution: sampling failed");
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BinomialDistribution<IntType>::variateWaiting(IntType number, G &g) const
{
    /// waiting algorithm, using
    /// sum of geometrically distributed variables
    IntType X = 0, sum = 0;
    while (true) {
        IntType wait = geometric(g);
        if (wait >= number - sum)
            return X;
        sum += wait + 1;
        ++X;
    }
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BinomialDistribution<IntType>::variateWaiting(IntType number, double probability, G &g)
{
    IntType X = 0, sum = 0;
    while (true) {
        IntType wait = GeometricRand<IntType>::Variate(probability, g);
        if (wait >= number - sum)
            return X;
        sum += wait + 1;
        ++X;
    }
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BinomialDistribution<IntType>::variateBernoulliSum(IntType number, double probability, G &g)
{
    IntType var = 0;
    if (math::detail::areClose(probability, 0.5)) {
        for (int i = 0; i != number; ++i)
            var += static_cast<IntType>(detail::randomBits(g, 1));
    } else {
        for (int i = 0; i != number; ++i)
            var += static_cast<IntType>(detail::uniform01<double>(g) < probability);
    }
    return var;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BinomialDistribution<IntType>::operator()(G &g) const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    switch (genId) {
    case WAITING: {
        IntType var = variateWaiting(n, g);
        return (p <= 0.5) ? var : n - var;
    }
    case REJECTION: {
        /// if X ~ Bin(n, p') and Y ~ Bin(n - X, (p - p') / (1 - p'))
        /// then Z = X + Y ~ Bin(n, p)
        IntType Z = variateRejection(g);
        if (pRes > 0)
            Z += variateWaiting(n - Z, g);
        return (p > 0.5) ? n - Z : Z;
    }
    case BERNOULLI_SUM:
        return variateBernoulliSum(n, p, g);
    default:
        throw std::invalid_argument("Binomial distribution: invalid generator id");
    }
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BinomialDistribution<IntType>::Variate(IntType number, double probability, G &g)
{
    /// sanity check
    if (!(number >= 0))
        throw std::invalid_argument("Binomial distribution: number should be positive, but it's equal to " +
                                    std::to_string(number));
    if (!(probability >= 0.0 && probability <= 1.0))
        throw std::invalid_argument(
            "Binomial distribution: probability parameter should in interval [0, 1], but it's equal to " +
            std::to_string(probability));
    if (probability == 0.0)
        return 0;
    if (probability == 1.0)
        return number;

    if (number < 10)
        return variateBernoulliSum(number, probability, g);
    if (probability < 0.5)
        return variateWaiting(number, probability, g);
    return number - variateWaiting(number, 1.0 - probability, g);
}

template <typename IntType>
void BinomialDistribution<IntType>::sample(detail::EngineRef g, std::span<IntType> outputData) const
{
    if (p == 0.0) {
        std::fill(outputData.begin(), outputData.end(), 0);
        return;
    }
    if (math::detail::areClose(p, 1.0)) {
        std::fill(outputData.begin(), outputData.end(), n);
        return;
    }

    GENERATOR_ID genId = GetIdOfUsedGenerator();
    switch (genId) {
    case WAITING: {
        if (p <= 0.5) {
            for (IntType &var : outputData)
                var = variateWaiting(n, g);
        } else {
            for (IntType &var : outputData)
                var = n - variateWaiting(n, g);
        }
        return;
    }
    case REJECTION: {
        for (IntType &var : outputData)
            var = variateRejection(g);
        if (pRes > 0) {
            for (IntType &var : outputData)
                var += variateWaiting(n - var, g);
        }
        if (p > 0.5) {
            for (IntType &var : outputData)
                var = n - var;
        }
        return;
    }
    case BERNOULLI_SUM:
    default: {
        for (IntType &var : outputData)
            var = variateBernoulliSum(n, p, g);
        return;
    }
    }
}

template <typename IntType>
long double BinomialDistribution<IntType>::Mean() const
{
    return np;
}

template <typename IntType>
long double BinomialDistribution<IntType>::Variance() const
{
    return np * q;
}

template <typename IntType>
IntType BinomialDistribution<IntType>::Median() const
{
    /// The median is floor(np) or ceil(np)
    IntType lower = std::floor(np);
    return (this->F(lower) >= 0.5) ? lower : lower + 1;
}

template <typename IntType>
IntType BinomialDistribution<IntType>::Mode() const
{
    return std::floor(np + p);
}

template <typename IntType>
long double BinomialDistribution<IntType>::Skewness() const
{
    return (q - p) / std::sqrt(np * q);
}

template <typename IntType>
long double BinomialDistribution<IntType>::ExcessKurtosis() const
{
    long double y = 1.0 / (p * q);
    y -= 6.0;
    return y / n;
}

template <typename IntType>
std::complex<double> BinomialDistribution<IntType>::CFImpl(double t) const
{
    std::complex<double> y(q + p * std::cos(t), p * std::sin(t));
    return std::pow(y, n);
}

template <typename IntType>
void BinomialDistribution<IntType>::FitProbability(std::span<const IntType> sample)
{
    if (!this->allElementsAreNonNegative(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->NON_NEGATIVITY_VIOLATION);
    if (!this->allElementsAreNotGreaterThan(n, sample))
        throw this->fitError(FitErrorCode::WrongSample, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(n));
    SetParameters(n, this->GetSampleMean(sample) / n);
}

template <typename IntType>
BetaRand<> BinomialDistribution<IntType>::FitProbabilityBayes(std::span<const IntType> sample,
                                                              const BetaDistribution<> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNonNegative(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->NON_NEGATIVITY_VIOLATION);
    if (!this->allElementsAreNotGreaterThan(n, sample))
        throw this->fitError(FitErrorCode::WrongSample, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(n));
    int N = sample.size();
    double sum = this->GetSampleSum(sample);
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    BetaRand posteriorDistribution(sum + alpha, N * n - sum + beta);
    SetParameters(n, MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template <typename IntType>
BetaRand<> BinomialDistribution<IntType>::FitProbabilityMinimax(std::span<const IntType> sample)
{
    double shape = 0.5 * std::sqrt(n);
    BetaRand B(shape, shape);
    return FitProbabilityBayes(sample, B);
}
template <typename IntType>
std::string BinomialRand<IntType>::Name() const
{
    return "Binomial(" + this->toStringWithPrecision(this->GetNumber()) + ", " +
           this->toStringWithPrecision(this->GetProbability()) + ")";
}

} // namespace randlib

RANDLIB_DETAIL_SUPPRESS_DEPRECATED_END

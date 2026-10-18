#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The HyperGeometricRand class <BR>
 * Hypergeometric distribution
 *
 * X ~ HG(N, K, n)
 */
template <typename IntType>
class HyperGeometricRand final : public DiscreteDistribution<IntType>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }

    IntType N = 1; ///< population size
    IntType K = 1; /// number of possible successes
    IntType n = 1; /// number of draws
    double pmfCoef = 1; ///< C(N, n)

public:
    HyperGeometricRand(IntType totalSize = 1, IntType drawsNum = 1, IntType successesNum = 1);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] IntType MinValue() const override { return std::max(static_cast<IntType>(0), n - N + K); }
    [[nodiscard]] IntType MaxValue() const override { return std::min(n, K); }

    void SetParameters(IntType totalSize, IntType drawsNum, IntType successesNum);
    [[nodiscard]] inline IntType GetTotalSize() { return N; }
    [[nodiscard]] inline IntType GetDrawsNum() { return n; }
    [[nodiscard]] inline IntType GetSuccessesNum() { return K; }

    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

private:
    /// below this many draws, after the symmetries, the urn is sampled directly
    static constexpr IntType SMALL_SAMPLE = 10;
    template <std::uniform_random_bit_generator G>
    IntType variateByUrn(IntType draws, IntType successes, G &g) const;
    template <std::uniform_random_bit_generator G>
    IntType variateRatioOfUniforms(IntType draws, IntType successes, G &g) const;

public:
    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
};

} // namespace randlib

#include <randlib/discrete/bernoulli.hpp>
#include <randlib/discrete/beta_binomial.hpp>

namespace randlib
{

template <typename IntType>
HyperGeometricRand<IntType>::HyperGeometricRand(IntType totalSize, IntType drawsNum, IntType successesNum)
{
    SetParameters(totalSize, drawsNum, successesNum);
}

template <typename IntType>
std::string HyperGeometricRand<IntType>::Name() const
{
    return "Hypergeometric(" + this->toStringWithPrecision(N) + ", " + this->toStringWithPrecision(n) + ", " +
           this->toStringWithPrecision(K) + ")";
}

template <typename IntType>
void HyperGeometricRand<IntType>::SetParameters(IntType totalSize, IntType drawsNum, IntType successesNum)
{
    if (!(totalSize > 0 && drawsNum > 0 && successesNum > 0))
        throw std::invalid_argument("HyperGeometric distribution: all parameters should be positive");
    if (!(drawsNum <= totalSize))
        throw std::invalid_argument("HyperGeometric distribution: total size should be greater than draws number");
    if (!(successesNum <= totalSize))
        throw std::invalid_argument("HyperGeometric distribution: total size should be greater than successes number");

    N = totalSize;
    n = drawsNum;
    K = successesNum;

    pmfCoef = math::lfact(K);
    pmfCoef += math::lfact(N - K);
    pmfCoef += math::lfact(N - n);
    pmfCoef += math::lfact(n);
    pmfCoef -= math::lfact(N);
}

template <typename IntType>
double HyperGeometricRand<IntType>::P(const IntType &k) const
{
    return (k < MinValue() || k > MaxValue()) ? 0.0 : std::exp(logP(k));
}

template <typename IntType>
double HyperGeometricRand<IntType>::logP(const IntType &k) const
{
    if (k < MinValue() || k > MaxValue())
        return -INFINITY;
    double y = math::lfact(k);
    y += math::lfact(K - k);
    y += math::lfact(n - k);
    y += math::lfact(N - K - n + k);
    return pmfCoef - y;
}

template <typename IntType>
double HyperGeometricRand<IntType>::F(const IntType &k) const
{
    IntType minVal = MinValue();
    if (k < minVal)
        return 0.0;
    IntType maxVal = MaxValue();
    if (k >= maxVal)
        return 1.0;
    if (k - minVal <= maxVal - k) {
        /// sum P(X = i) going forward until k
        double sum = 0;
        for (IntType i = minVal; i <= k; ++i)
            sum += P(i);
        return sum;
    }
    /// going backwards is faster
    double sum = 1.0;
    for (IntType i = k + 1; i <= maxVal; ++i)
        sum -= P(i);
    return sum;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType HyperGeometricRand<IntType>::operator()(G &g) const
{
    /// X(N, n, K) has the law of X(N, K, n), K - X(N, N - n, K) and n - X(N, n, N - K): the draws and the
    /// successes are both reduced to at most half of the population first
    const bool complementDraws = 2 * n > N, complementSuccesses = 2 * K > N;
    const IntType draws = complementDraws ? N - n : n;
    const IntType successes = complementSuccesses ? N - K : K;
    const IntType reduced =
        (draws < SMALL_SAMPLE) ? variateByUrn(draws, successes, g) : variateRatioOfUniforms(draws, successes, g);
    IntType result = complementSuccesses ? draws - reduced : reduced;
    if (complementDraws)
        result = K - result;
    return result;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType HyperGeometricRand<IntType>::variateByUrn(IntType draws, IntType successes, G &g) const
{
    /// sequential draws without replacement, at most draws < SMALL_SAMPLE of them
    IntType sum = 0;
    for (IntType i = 0; i < draws && sum < successes; ++i) {
        const double p = static_cast<double>(successes - sum) / (N - i);
        if (BernoulliRand::Variate(p, g))
            ++sum;
    }
    return sum;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType HyperGeometricRand<IntType>::variateRatioOfUniforms(IntType draws, IntType successes, G &g) const
{
    /// HRUA: ratio of uniforms under Stadlober's table-mountain hat (Stadlober, 1989 and 1990), whose expected
    /// number of iterations is bounded for all parameters; draws, successes ≤ N / 2
    static constexpr double D1 = 1.7155277699214135, D2 = 0.8989161620588988;
    const IntType failures = N - successes;
    const double p = static_cast<double>(successes) / N, q = static_cast<double>(failures) / N;
    const double a = draws * p + 0.5;
    const double c = std::sqrt(static_cast<double>(N - draws) * draws * p * q / (N - 1) + 0.5);
    const double h = D1 * c + D2;
    const IntType mode = static_cast<IntType>(std::floor((draws + 1.0) * (successes + 1.0) / (N + 2.0)));
    /// log(P(k) / P(mode)) as a sum of log-factorial ratios, free of the cancellation of large log-factorials
    auto logRatioToPeak = [draws, successes, failures, mode](IntType k) {
        auto ratio = [](IntType x, IntType y) {
            return math::logFactorialRatio(static_cast<size_t>(x), static_cast<size_t>(y));
        };
        return ratio(mode, k) + ratio(successes - mode, successes - k) + ratio(draws - mode, draws - k) +
               ratio(failures - draws + mode, failures - draws + k);
    };
    const double end = static_cast<double>(std::min(draws, successes)) + 1;
    for (size_t iter = 0; iter < 1000 * ProbabilityDistribution<IntType>::MAX_ITER_REJECTION; ++iter) {
        const double u = detail::uniform01<double>(g), v = detail::uniform01<double>(g);
        const double x = a + h * (v - 0.5) / u;
        if (x < 0 || x >= end)
            continue;
        const IntType k = static_cast<IntType>(std::floor(x));
        const double t = logRatioToPeak(k);
        if (u * (4 - u) - 3 <= t)
            return k;
        if (u * (u - t) >= 1)
            continue;
        if (2 * std::log(u) <= t)
            return k;
    }
    throw std::runtime_error("HyperGeometric distribution: sampling failed");
}

template <typename IntType>
long double HyperGeometricRand<IntType>::Mean() const
{
    return static_cast<long double>(n) * K / N;
}

template <typename IntType>
long double HyperGeometricRand<IntType>::Variance() const
{
    if (n == N || K == N)
        return 0.0L;
    long double numerator = n;
    numerator *= K;
    numerator *= N - K;
    numerator *= N - n;
    long double denominator = N;
    denominator *= N;
    denominator *= N - 1;
    return numerator / denominator;
}

template <typename IntType>
IntType HyperGeometricRand<IntType>::Mode() const
{
    long double mode = (static_cast<long double>(n) + 1.0L) * (static_cast<long double>(K) + 1.0L);
    return static_cast<IntType>(std::floor(mode / (static_cast<long double>(N) + 2.0L)));
}

template <typename IntType>
long double HyperGeometricRand<IntType>::Skewness() const
{
    if (n == N || K == N)
        return NAN;
    if (N == 2)
        return 0.0L;
    long double population = N, draws = n, successes = K;
    long double skewness = population - 1.0L;
    skewness /= draws;
    skewness /= successes;
    skewness /= population - successes;
    skewness /= population - draws;
    skewness = std::sqrt(skewness);
    skewness *= population - 2.0L * successes;
    skewness *= population - 2.0L * draws;
    return skewness / (population - 2.0L);
}

template <typename IntType>
long double HyperGeometricRand<IntType>::ExcessKurtosis() const
{
    long double variance = Variance();
    if (variance == 0.0L)
        return NAN;
    if (N <= 3) {
        long double mean = Mean(), fourthMoment = 0.0L;
        for (IntType k = MinValue(); k <= MaxValue(); ++k) {
            long double centered = k - mean;
            fourthMoment += centered * centered * centered * centered * P(k);
        }
        return fourthMoment / (variance * variance) - 3.0L;
    }

    long double population = N, draws = n, successes = K;
    long double numerator = population * (population + 1.0L);
    long double a1 = successes * (population - successes);
    long double a2 = draws * (population - draws);
    numerator -= 6 * (a1 + a2);
    numerator *= population * population * (population - 1.0L);
    long double aux = draws * successes * (population - successes) * (population - draws);
    numerator += 6.0L * aux * (5.0L * population - 6.0L);
    long double denominator = aux * (population - 2.0L) * (population - 3.0L);
    return numerator / denominator;
}

} // namespace randlib

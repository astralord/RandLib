#pragma once

#include <functional>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/beta.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The BetaBinomialRand class <BR>
 * Beta-Binomial distribution
 *
 * Notation: X ~ BB(n, α, β)
 *
 * Related distributions: <BR>
 * If X ~ Binomial(n, p), where p ~ Beta(α, β), then X ~ BB(n, α, β)
 */
template <typename IntType>
class BetaBinomialRand final : public DiscreteDistribution<IntType>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<IntType> outputData) const override;

    IntType n = 1; ///< number of experiments
    double pmfCoef = 0; ///< log(n!) - log(Γ(α + β + n)) - log(B(α, β))
    BetaRand<double> B{};

public:
    BetaBinomialRand(IntType number = 1, double shape1 = 1, double shape2 = 1);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] IntType MinValue() const override { return 0; }
    [[nodiscard]] IntType MaxValue() const override { return n; }

    void SetParameters(IntType number, double shape1, double shape2);
    [[nodiscard]] inline IntType GetNumber() const noexcept { return n; }
    [[nodiscard]] inline double GetAlpha() const { return B.GetAlpha(); }
    [[nodiscard]] inline double GetBeta() const { return B.GetBeta(); }

    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;

private:
    template <std::uniform_random_bit_generator G>
    IntType VariateUniform(G &g) const;
    template <std::uniform_random_bit_generator G>
    IntType VariateBeta(G &g) const;

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
};

} // namespace randlib

#include <randlib/discrete/binomial.hpp>
#include <randlib/discrete/uniform_discrete.hpp>

namespace randlib
{

template <typename IntType>
BetaBinomialRand<IntType>::BetaBinomialRand(IntType number, double shape1, double shape2)
{
    SetParameters(number, shape1, shape2);
}

template <typename IntType>
std::string BetaBinomialRand<IntType>::Name() const
{
    return "Beta-Binomial(" + this->toStringWithPrecision(GetNumber()) + ", " +
           this->toStringWithPrecision(GetAlpha()) + ", " + this->toStringWithPrecision(GetBeta()) + ")";
}

template <typename IntType>
void BetaBinomialRand<IntType>::SetParameters(IntType number, double shape1, double shape2)
{
    if (!(shape1 > 0.0 && shape2 > 0.0))
        throw std::invalid_argument("Beta-Binomial distribution: shape parameters should be positive");
    if (!(number > 0))
        throw std::invalid_argument("Beta-Binomial distribution: number should be positive");
    n = number;
    B.SetShapes(shape1, shape2);
    pmfCoef = math::lfact(n);
    pmfCoef -= std::lgammal(B.GetAlpha() + B.GetBeta() + n);
    pmfCoef -= B.GetLogBetaFunction();
}

template <typename IntType>
double BetaBinomialRand<IntType>::P(const IntType &k) const
{
    return (k < 0 || k > n) ? 0.0 : std::exp(logP(k));
}

template <typename IntType>
double BetaBinomialRand<IntType>::logP(const IntType &k) const
{
    if (k < 0 || k > n)
        return -INFINITY;
    double y = std::lgammal(k + B.GetAlpha());
    y += std::lgammal(n - k + B.GetBeta());
    y -= math::lfact(k);
    y -= math::lfact(n - k);
    return pmfCoef + y;
}

template <typename IntType>
double BetaBinomialRand<IntType>::F(const IntType &k) const
{
    if (k < 0)
        return 0.0;
    if (k >= n)
        return 1.0;
    double sum = 0.0;
    int i = 0;
    do {
        sum += P(i);
    } while (++i <= k);
    return sum;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BetaBinomialRand<IntType>::VariateUniform(G &g) const
{
    return UniformDiscreteRand<IntType>::StandardVariate(0, n, g);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BetaBinomialRand<IntType>::VariateBeta(G &g) const
{
    double p = B(g);
    return BinomialDistribution<IntType>::Variate(n, p, g);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType BetaBinomialRand<IntType>::operator()(G &g) const
{
    return (B.GetAlpha() == 1 && B.GetBeta() == 1) ? VariateUniform(g) : VariateBeta(g);
}

template <typename IntType>
void BetaBinomialRand<IntType>::sample(detail::EngineRef g, std::span<IntType> outputData) const
{
    if (B.GetAlpha() == 1 && B.GetBeta() == 1) {
        for (IntType &var : outputData)
            var = VariateUniform(g);
    } else {
        for (IntType &var : outputData)
            var = VariateBeta(g);
    }
}

template <typename IntType>
long double BetaBinomialRand<IntType>::Mean() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    return n * alpha / (alpha + beta);
}

template <typename IntType>
long double BetaBinomialRand<IntType>::Variance() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    double alphaPBeta = alpha + beta;
    double numerator = n * alpha * beta * (alphaPBeta + n);
    double denominator = alphaPBeta * alphaPBeta;
    denominator *= (alphaPBeta + 1);
    return numerator / denominator;
}

template <typename IntType>
IntType BetaBinomialRand<IntType>::Mode() const
{
    IntType mode = (IntType)(n * B.Mode());
    double logPmode = this->logP(mode);
    if (this->logP(mode + 1) > logPmode)
        return mode + 1;
    if (this->logP(mode - 1) > logPmode)
        return mode - 1;
    return mode;
}

template <typename IntType>
long double BetaBinomialRand<IntType>::Skewness() const
{
    long double alpha = B.GetAlpha();
    long double beta = B.GetBeta();
    long double alphaPBeta = alpha + beta;
    long double res = (1 + alphaPBeta) / (n * alpha * beta * (alphaPBeta + n));
    res = std::sqrt(res);
    res *= (alphaPBeta + 2 * n) * (beta - alpha);
    res /= alphaPBeta + 2;
    return res;
}

template <typename IntType>
long double BetaBinomialRand<IntType>::ExcessKurtosis() const
{
    long double alpha = B.GetAlpha();
    long double beta = B.GetBeta();
    long double alphaPBeta = alpha + beta;
    long double alphaBetaN = alpha * beta * n;
    long double res = alpha * beta * (n - 2);
    res += 2 * (double)n * n;
    res -= alphaBetaN * (6 - n) / alphaPBeta;
    res -= 6 * alphaBetaN * n / (alphaPBeta * alphaPBeta);
    res *= 3;
    res += alphaPBeta * (alphaPBeta - 1 + 6 * n);
    res *= alphaPBeta * alphaPBeta * (1 + alphaPBeta);
    res /= (alphaBetaN * (alphaPBeta + 2) * (alphaPBeta + 3) * (alphaPBeta + n));
    return res - 3.0;
}

} // namespace randlib

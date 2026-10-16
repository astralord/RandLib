#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The UniformDiscreteRand class <BR>
 * Uniform discrete distribution
 *
 * Notation: X ~ U{a, ..., b}
 *
 * P(X = k) = 1 / (b - a + 1) for a <= k <= b
 */
template <typename IntType>
class UniformDiscreteRand final : public DiscreteDistribution<IntType>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }

    size_t n = 1; ///< number of possible outcomes
    IntType a = 0; ///< min bound
    IntType b = 0; ///< max bound
    double nInv = 1; ///< 1/n
    double logN = 0; ///< log(n)

public:
    UniformDiscreteRand(IntType minValue = 0, IntType maxValue = 1);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] IntType MinValue() const override { return a; }
    [[nodiscard]] IntType MaxValue() const override { return b; }

    void SetBoundaries(IntType minValue, IntType maxValue);

    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static IntType StandardVariate(IntType minValue, IntType maxValue, G &g);

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Median() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Entropy() const;
    [[nodiscard]] double LikelihoodFunction(std::span<const IntType> sample) const override;
    [[nodiscard]] double LogLikelihoodFunction(std::span<const IntType> sample) const override;

    /**
     * fit bounds via maximum-likelihood method
     * @param sample
     */
    RANDLIB_DEPRECATED_FIT_MEMBER void Fit(std::span<const IntType> sample);
};

} // namespace randlib

namespace randlib
{

template <typename IntType>
UniformDiscreteRand<IntType>::UniformDiscreteRand(IntType minValue, IntType maxValue)
{
    SetBoundaries(minValue, maxValue);
}

template <typename IntType>
std::string UniformDiscreteRand<IntType>::Name() const
{
    return "Uniform Discrete(" + this->toStringWithPrecision(MinValue()) + ", " +
           this->toStringWithPrecision(MaxValue()) + ")";
}

template <typename IntType>
void UniformDiscreteRand<IntType>::SetBoundaries(IntType minValue, IntType maxValue)
{
    if (!(minValue <= maxValue))
        throw std::invalid_argument(
            "Uniform discrete distribution: minimum value shouldn't be greater than maximum value");

    a = minValue;
    b = maxValue;

    n = b - a + 1;
    nInv = 1.0 / n;
    logN = std::log(n);
}

template <typename IntType>
double UniformDiscreteRand<IntType>::P(const IntType &k) const
{
    return (k < a || k > b) ? 0.0 : nInv;
}

template <typename IntType>
double UniformDiscreteRand<IntType>::logP(const IntType &k) const
{
    return (k < a || k > b) ? -INFINITY : -logN;
}

template <typename IntType>
double UniformDiscreteRand<IntType>::F(const IntType &k) const
{
    if (k < a)
        return 0.0;
    if (k > b)
        return 1.0;
    return (k - a + 1) * nInv;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType UniformDiscreteRand<IntType>::operator()(G &g) const
{
    return a + static_cast<IntType>(detail::uniformBelow(g, n));
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType UniformDiscreteRand<IntType>::StandardVariate(IntType minValue, IntType maxValue, G &g)
{
    if (maxValue <= minValue)
        return minValue;
    auto n = static_cast<std::uint64_t>(maxValue - minValue) + 1;
    return minValue + static_cast<IntType>(detail::uniformBelow(g, n));
}

template <typename IntType>
long double UniformDiscreteRand<IntType>::Mean() const
{
    return 0.5 * (b + a);
}

template <typename IntType>
long double UniformDiscreteRand<IntType>::Variance() const
{
    double nm1 = n - 1;
    double np1 = n + 1;
    return nm1 * np1 / 12;
}

template <typename IntType>
IntType UniformDiscreteRand<IntType>::Median() const
{
    return (b + a) >> 1;
}

template <typename IntType>
IntType UniformDiscreteRand<IntType>::Mode() const
{
    /// this can be any value in [a, b]
    return 0.5 * (a + b);
}

template <typename IntType>
long double UniformDiscreteRand<IntType>::Skewness() const
{
    return 0.0;
}

template <typename IntType>
long double UniformDiscreteRand<IntType>::ExcessKurtosis() const
{
    double kurt = n;
    kurt *= n;
    --kurt;
    kurt = 2.0 / kurt;
    ++kurt;
    return -1.2 * kurt;
}

template <typename IntType>
std::complex<double> UniformDiscreteRand<IntType>::CFImpl(double t) const
{
    double at = a * t;
    double bp1t = (b + 1) * t;
    double reNum = std::cos(at) - std::cos(bp1t);
    double imNum = std::sin(at) - std::sin(bp1t);
    std::complex<double> numerator(reNum, imNum);
    std::complex<double> denominator(1.0 - std::cos(t), -std::sin(t));
    return nInv * numerator / denominator;
}

template <typename IntType>
long double UniformDiscreteRand<IntType>::Entropy() const
{
    return logN;
}

template <typename IntType>
double UniformDiscreteRand<IntType>::LikelihoodFunction(std::span<const IntType> sample) const
{
    bool sampleIsInsideInterval =
        this->allElementsAreNotSmallerThan(a, sample) && this->allElementsAreNotGreaterThan(b, sample);
    return sampleIsInsideInterval ? std::pow(n, -static_cast<double>(sample.size())) : 0.0;
}

template <typename IntType>
double UniformDiscreteRand<IntType>::LogLikelihoodFunction(std::span<const IntType> sample) const
{
    bool sampleIsInsideInterval =
        this->allElementsAreNotSmallerThan(a, sample) && this->allElementsAreNotGreaterThan(b, sample);
    int sample_size = sample.size();
    return sampleIsInsideInterval ? -sample_size * logN : -INFINITY;
}

template <typename IntType>
void UniformDiscreteRand<IntType>::Fit(std::span<const IntType> sample)
{
    IntType minVar = *std::min_element(sample.begin(), sample.end());
    IntType maxVar = *std::max_element(sample.begin(), sample.end());
    this->SetBoundaries(minVar, maxVar);
}

} // namespace randlib

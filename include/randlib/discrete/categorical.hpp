#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The CategoricalRand class <BR>
 *
 * P(X = k) = p_k for k = {0, ..., K-1}
 *
 * Notation: X ~ Cat(p_0, ..., p_{K-1})
 *
 * Related distributions: <BR>
 * X ~ Multin(1, p_0, ..., p_{K-1}) <BR>
 * If X ~ Bernoulli(p), then X ~ Cat(1 - p, p) <BR>
 * If X ~ Uniform-Discrete(0, K), then X ~ Cat(p, ..., p) with p = 1 / (K + 1)
 */
template <typename IntType>
class CategoricalRand final : public DiscreteDistribution<IntType>
{
    std::vector<double> prob{1.0}; ///< vector of probabilities
    IntType K = 1; ///< number of possible outcomes

public:
    explicit CategoricalRand(std::vector<double>&& probabilities = {1.0});
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] IntType MinValue() const override { return 0; }
    [[nodiscard]] IntType MaxValue() const override { return K - 1; }

    void SetProbabilities(std::vector<double>&& probabilities);
    [[nodiscard]] std::vector<double> GetProbabilities() { return prob; }

    [[nodiscard]] double P(const IntType & k) const override;
    [[nodiscard]] double logP(const IntType & k) const override;
    [[nodiscard]] double F(const IntType & k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;
protected:
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Mode() const override;

private:
    IntType quantileImpl(double p) const override;
    IntType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template < typename IntType >
CategoricalRand<IntType>::CategoricalRand(std::vector<double>&& probabilities)
{
    SetProbabilities(std::move(probabilities));
}

template < typename IntType >
std::string CategoricalRand<IntType>::Name() const
{
    std::string str = "Categorical(";
    for (int i = 0; i != K - 1; ++i)
        str += this->toStringWithPrecision(prob[i]) + ", ";
    return str + this->toStringWithPrecision(prob[K - 1]) + ")";
}

template < typename IntType >
void CategoricalRand<IntType>::SetProbabilities(std::vector<double> &&probabilities)
{
    if (probabilities.empty())
        throw std::invalid_argument("Categorical distribution: there should be at least one probability");
    for (double probability : probabilities) {
        if (!(probability >= 0.0))
            throw std::invalid_argument("Categorical distribution: probabilities should be non-negative");
    }
    double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    if (!(std::fabs(sum - 1.0) <= 1e-12 * probabilities.size()))
        throw std::invalid_argument("Categorical distribution: probability parameters should sum to 1");
    prob = std::move(probabilities);
    K = prob.size();
}

template < typename IntType >
double CategoricalRand<IntType>::P(const IntType & k) const
{
    return (k < 0 || k >= K) ? 0.0 : prob[k];
}

template < typename IntType >
double CategoricalRand<IntType>::logP(const IntType & k) const
{
    return std::log(P(k));
}

template < typename IntType >
double CategoricalRand<IntType>::F(const IntType & k) const
{
    if (k < 0)
        return 0.0;
    if (k >= K)
        return 1.0;
    if (2 * k <= K) {
        double sum = 0.0;
        for (int i = 0; i <= k; ++i)
            sum += prob[i];
        return sum;
    }
    double sum = 1.0;
    for (int i = K - 1; i > k; --i)
        sum -= prob[i];
    return sum;
}

template < typename IntType >
template <std::uniform_random_bit_generator G>
IntType CategoricalRand<IntType>::operator()(G &g) const
{
    double U = detail::uniform01<double>(g);
    return quantileImpl(U);
}

template < typename IntType >
long double CategoricalRand<IntType>::Mean() const
{
    long double sum = 0.0;
    for (int i = 1; i != K; ++i)
        sum += i * prob[i];
    return sum;
}

template < typename IntType >
long double CategoricalRand<IntType>::Variance() const
{
    long double mean = 0.0, secMom = 0.0;
    for (int i = 1; i != K; ++i) {
        long double aux = i * prob[i];
        mean += aux;
        secMom += i * aux;
    }
    return secMom - mean * mean;
}

template < typename IntType >
IntType CategoricalRand<IntType>::Mode() const
{
    auto maxProbIt = std::max_element(prob.begin(), prob.end());
    return std::distance(prob.begin(), maxProbIt);
}

template < typename IntType >
IntType CategoricalRand<IntType>::quantileImpl(double p) const
{
    /// smallest k with F(k) ≥ p
    double sum = 0.0;
    for (IntType i = 0; i != K; ++i) {
        sum += prob[i];
        if (sum >= p)
            return i;
    }
    return K - 1;
}

template < typename IntType >
IntType CategoricalRand<IntType>::quantileImpl1m(double p) const
{
    /// smallest k with S(k) ≤ p, where S(k - 1) = S(k) + P(k)
    double tail = 0.0;
    for (IntType k = K - 1; k > 0; --k) {
        tail += prob[k];
        if (tail > p)
            return k;
    }
    return 0;
}

template < typename IntType >
std::complex<double> CategoricalRand<IntType>::CFImpl(double t) const
{
    double re = 0.0;
    double im = 0.0;
    for (int i = 0; i != K; ++i) {
        re += prob[i] * std::cos(t * i);
        im += prob[i] * std::sin(t * i);
    }
    return std::complex<double>(re, im);
}

} // namespace randlib

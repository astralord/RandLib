#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/discrete/binomial.hpp>

namespace randlib
{

/**
 * @brief The BernoulliRand class <BR>
 * Bernoulli distribution
 *
 * P(X = k) = p * 1_{k = 1} + (1 - p) * 1_{k = 0}
 *
 * Notation: X ~ Bernoulli(p)
 *
 * Related distributions: <BR>
 * X ~ Binomial(1, p) <BR>
 * X ~ Multin(1, 1 - p, p) <BR>
 * 2X - 1 ~ Rademacher
 */
class BernoulliRand final : public BinomialDistribution<int>
{
    int variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<int> outputData) const override;

public:
    explicit BernoulliRand(double probability = 0.5);
    [[nodiscard]] std::string Name() const override;

    void SetProbability(double probability);

    [[nodiscard]] double P(const int &k) const override;
    [[nodiscard]] double logP(const int &k) const override;
    [[nodiscard]] double F(const int &k) const override;
    [[nodiscard]] double S(const int &k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] int operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static int Variate(double probability, G &g);
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static int StandardVariate(G &g);

    [[nodiscard]] inline double Entropy();
};

} // namespace randlib

namespace randlib
{

inline BernoulliRand::BernoulliRand(double probability) : BinomialDistribution(1, probability)
{
}

inline std::string BernoulliRand::Name() const
{
    return "Bernoulli(" + this->toStringWithPrecision(GetProbability()) + ")";
}

inline void BernoulliRand::SetProbability(double probability)
{
    if (!(probability >= 0.0 && probability <= 1.0))
        throw std::invalid_argument("Bernoulli distribution: probability parameter should in interval [0, 1]");
    SetParameters(1, probability);
}

inline double BernoulliRand::P(const int &k) const
{
    return (k == 0) ? q : ((k == 1) ? p : 0);
}

inline double BernoulliRand::logP(const int &k) const
{
    return (k == 0) ? log1mProb : ((k == 1) ? logProb : -INFINITY);
}

inline double BernoulliRand::F(const int &k) const
{
    return (k < 0) ? 0.0 : ((k < 1) ? q : 1);
}

inline double BernoulliRand::S(const int &k) const
{
    return (k < 0) ? 1.0 : ((k < 1) ? p : 0.0);
}

template <std::uniform_random_bit_generator G>
int BernoulliRand::operator()(G &g) const
{
    return detail::uniform01<double>(g) < p;
}

template <std::uniform_random_bit_generator G>
int BernoulliRand::Variate(double probability, G &g)
{
    if (!(probability >= 0.0 && probability <= 1.0))
        throw std::invalid_argument("Bernoulli distribution: probability parameter should be in interval [0, 1]");
    return detail::uniform01<double>(g) < probability;
}

template <std::uniform_random_bit_generator G>
int BernoulliRand::StandardVariate(G &g)
{
    return static_cast<int>(detail::randomBits(g, 1));
}

inline void BernoulliRand::sample(detail::EngineRef g, std::span<int> outputData) const
{
    if (p == 0.5) {
        for (int &var : outputData)
            var = StandardVariate(g);
    } else {
        for (int &var : outputData)
            var = (*this)(g);
    }
}

inline double BernoulliRand::Entropy()
{
    return -(p * logProb + q * log1mProb);
}

} // namespace randlib

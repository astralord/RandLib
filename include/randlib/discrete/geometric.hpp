#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/discrete/negative_binomial.hpp>

namespace randlib
{

/**
 * @brief The GeometricRand class <BR>
 * Geometric distribution
 *
 * P(X = k) = p (1 - p)^k
 *
 * Notation: X ~ Geometric(p)
 *
 * Related distributions: <BR>
 * X ~ NB(1, p)
 */
template <typename IntType>
class GeometricRand final : public PascalRand<IntType>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<IntType> outputData) const override;

public:
    explicit GeometricRand(double probability = 0.5) : PascalRand<IntType>(1, probability) {}
    [[nodiscard]] std::string Name() const override;

    void SetProbability(double probability);

    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;
    [[nodiscard]] double S(const IntType &k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static IntType Variate(double probability, G &g);

    [[nodiscard]] IntType Median() const override;

    [[nodiscard]] long double Entropy() const;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/exponential.hpp>

namespace randlib
{

template <typename IntType>
std::string GeometricRand<IntType>::Name() const
{
    return "Geometric(" + this->toStringWithPrecision(this->GetProbability()) + ")";
}

template <typename IntType>
void GeometricRand<IntType>::SetProbability(double probability)
{
    if (!(probability >= 0.0 && probability <= 1.0))
        throw std::invalid_argument("Geometric distribution: probability parameter should be in interval [0, 1]");
    this->SetParameters(1, probability);
}

template <typename IntType>
double GeometricRand<IntType>::P(const IntType &k) const
{
    return (k < 0) ? 0 : this->p * std::exp(k * this->log1mProb);
}

template <typename IntType>
double GeometricRand<IntType>::logP(const IntType &k) const
{
    return (k < 0) ? -INFINITY : this->logProb + k * this->log1mProb;
}

template <typename IntType>
double GeometricRand<IntType>::F(const IntType &k) const
{
    return (k < 0) ? 0 : -std::expm1l((k + 1) * this->log1mProb);
}

template <typename IntType>
double GeometricRand<IntType>::S(const IntType &k) const
{
    return (k < 0) ? 1 : std::exp((k + 1) * this->log1mProb);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType GeometricRand<IntType>::operator()(G &g) const
{
    typename PascalRand<IntType>::GENERATOR_ID genId = this->GetIdOfUsedGenerator();
    if (genId == this->EXPONENTIAL)
        return this->variateGeometricThroughExponential(g);
    if (genId == this->TABLE)
        return this->variateGeometricByTable(g);
    /// unexpected return
    throw std::runtime_error("Geometric distribution: sampling failed");
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType GeometricRand<IntType>::Variate(double probability, G &g)
{
    if (!(probability <= 1.0 && probability >= 0.0))
        throw std::invalid_argument(
            "Geometric distribution: probability parameter should be in interval [0, 1], but it's equal to " +
            std::to_string(probability));

    /// here we use 0.05 instead of 0.08 because log(q) wasn't hashed
    if (probability < 0.05) {
        double rate = -std::log1pl(-probability);
        float X = detail::standardExponential<float>(g) / rate;
        return std::floor(X);
    }

    double U = detail::uniform01<double>(g);
    int x = 0;
    double prod = probability, sum = prod, qprob = 1.0 - probability;
    while (U > sum) {
        prod *= qprob;
        sum += prod;
        ++x;
    }
    return x;
}

template <typename IntType>
void GeometricRand<IntType>::sample(detail::EngineRef g, std::span<IntType> outputData) const
{
    typename PascalRand<IntType>::GENERATOR_ID genId = this->GetIdOfUsedGenerator();
    if (genId == this->EXPONENTIAL) {
        for (IntType &var : outputData)
            var = this->variateGeometricThroughExponential(g);
    } else if (genId == this->TABLE) {
        for (IntType &var : outputData)
            var = this->variateGeometricByTable(g);
    }
}

template <typename IntType>
IntType GeometricRand<IntType>::Median() const
{
    double median = -randlib::consts::ln2 / this->log1mProb;
    double flooredMedian = std::floor(median);
    return (math::detail::areClose(median, flooredMedian, 1e-8)) ? flooredMedian - 1 : flooredMedian;
}

template <typename IntType>
long double GeometricRand<IntType>::Entropy() const
{
    double a = -this->q * this->log1mProb;
    double b = -this->p * this->logProb;
    return (a + b) / (randlib::consts::ln2 * this->p);
}

} // namespace randlib

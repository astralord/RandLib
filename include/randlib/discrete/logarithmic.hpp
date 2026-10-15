#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The LogarithmicRand class <BR>
 * Logarithmic distribution
 *
 * P(X = k) = -p^k / [k log(1 - p)]
 *
 * X ~ Log(p)
 */
template <typename IntType>
class LogarithmicRand final : public DiscreteDistribution<IntType>
{
    double p = 0.5; ///< parameter of distribution
    double logProb = -randlib::consts::ln2; ///< log(p)
    double log1mProb = -randlib::consts::ln2; ///< log(q)
public:
    explicit LogarithmicRand(double probability);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] IntType MinValue() const override { return 1; }
    [[nodiscard]] IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetProbability(double probability);
    [[nodiscard]] inline double GetProbability() const noexcept { return p; }

    [[nodiscard]] double P(const IntType &k) const override;
    [[nodiscard]] double logP(const IntType &k) const override;

private:
    /**
     * @fn betaFun
     * @param a
     * @return B(p, a, 0), where B(x, a, b) denotes incomplete beta function,
     * using series expansion (converges for x < 1)
     */
    double betaFun(IntType a) const;

public:
    [[nodiscard]] double F(const IntType &k) const override;
    [[nodiscard]] double S(const IntType &k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

protected:
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }

public:
    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template <typename IntType>
LogarithmicRand<IntType>::LogarithmicRand(double probability)
{
    SetProbability(probability);
}

template <typename IntType>
std::string LogarithmicRand<IntType>::Name() const
{
    return "Logarithmic(" + this->toStringWithPrecision(GetProbability()) + ")";
}

template <typename IntType>
void LogarithmicRand<IntType>::SetProbability(double probability)
{
    if (!(probability > 0.0 && probability < 1.0))
        throw std::invalid_argument("Logarithmic distribution: probability parameter should in interval (0, 1)");
    p = probability;
    logProb = std::log(p);
    log1mProb = std::log1pl(-p);
}

template <typename IntType>
double LogarithmicRand<IntType>::P(const IntType &k) const
{
    return (k < 1) ? 0.0 : -std::pow(p, k) / (k * log1mProb);
}

template <typename IntType>
double LogarithmicRand<IntType>::logP(const IntType &k) const
{
    return (k < 1) ? -INFINITY : k * logProb - std::log(-k * log1mProb);
}

template <typename IntType>
double LogarithmicRand<IntType>::betaFun(IntType a) const
{
    /// Σ_{j≥0} p^j / (a + j), the first three terms written out
    double denom = a + 2;
    double sum = p * p / (a + 2) + p / (a + 1) + 1.0 / a;
    double add = 1;
    int i = 3;
    do {
        add = std::exp(i * logProb) / (++denom);
        sum += add;
        ++i;
    } while (add > MIN_POSITIVE * sum);
    return std::exp(a * logProb) * sum;
}

template <typename IntType>
double LogarithmicRand<IntType>::F(const IntType &k) const
{
    return (k < 1) ? 0.0 : 1 + betaFun(k + 1) / log1mProb;
}

template <typename IntType>
double LogarithmicRand<IntType>::S(const IntType &k) const
{
    return (k < 1) ? 1.0 : -betaFun(k + 1) / log1mProb;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType LogarithmicRand<IntType>::operator()(G &g) const
{
    /// Kemp's second accelerated generator
    /// p. 548, "Non-Uniform Random Variate Generation" by Luc Devroye
    float V = detail::uniform01<float>(g);
    if (V >= p)
        return 1.0;
    float U = detail::uniform01<float>(g);
    double y = -std::expm1l(U * log1mProb);
    if (V > y)
        return 1.0;
    if (V > y * y)
        return 2.0;
    return std::floor(1.0 + std::log(V) / std::log(y));
}

template <typename IntType>
long double LogarithmicRand<IntType>::Mean() const
{
    return -p / (1.0 - p) / log1mProb;
}

template <typename IntType>
long double LogarithmicRand<IntType>::Variance() const
{
    long double var = p / log1mProb + 1;
    var /= log1mProb;
    var *= p;
    long double q = 1.0 - p;
    return -var / (q * q);
}

template <typename IntType>
std::complex<double> LogarithmicRand<IntType>::CFImpl(double t) const
{
    std::complex<double> y(std::cos(t), std::sin(t));
    y = std::log(1.0 - p * y);
    return y / log1mProb;
}

template <typename IntType>
IntType LogarithmicRand<IntType>::Mode() const
{
    return 1.0;
}

} // namespace randlib

#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>
#include <randlib/continuous/pareto.hpp>

namespace randlib
{

/**
 * @brief The YuleRand class <BR>
 * Yule distribution
 *
 * Notation: X ~ Yule(ρ)
 *
 * Related distributions: <BR>
 * If Y ~ Pareto(ρ, 1) and Z ~ Geometric(1 / Y), then Z + 1 ~ Yule(ρ)
 */
template <typename IntType>
class YuleRand final : public DiscreteDistribution<IntType>
{
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }

    double rho = 0; ///< shape ρ
    double lgamma1pRo = 0; /// log(Γ(1 + ρ))

    ParetoRand<double> X;

public:
    explicit YuleRand(double shape);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] IntType MinValue() const override { return 1; }
    [[nodiscard]] IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetShape(double shape);
    [[nodiscard]] inline double GetShape() const noexcept { return rho; }

    [[nodiscard]] double logP(const IntType &k) const override;
    [[nodiscard]] double F(const IntType &k) const override;
    [[nodiscard]] double S(const IntType &k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;

    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static IntType Variate(double shape, G &g);

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    template <std::uniform_random_bit_generator G>
    static IntType variateGivenPareto(double pareto, G &g);
};

} // namespace randlib

#include <randlib/discrete/geometric.hpp>
#include <randlib/continuous/exponential.hpp>

namespace randlib
{

template <typename IntType>
YuleRand<IntType>::YuleRand(double shape) : X(shape, 1.0)
{
    SetShape(shape);
}

template <typename IntType>
std::string YuleRand<IntType>::Name() const
{
    return "Yule(" + this->toStringWithPrecision(GetShape()) + ")";
}

template <typename IntType>
void YuleRand<IntType>::SetShape(double shape)
{
    if (!(shape > 0.0))
        throw std::invalid_argument("Yule distribution: shape should be positive");
    rho = shape;
    lgamma1pRo = std::lgammal(rho + 1);
    X.SetShape(rho);
}

template <typename IntType>
double YuleRand<IntType>::logP(const IntType &k) const
{
    if (k < 1)
        return -INFINITY;
    double y = lgamma1pRo;
    y += math::lfact(k - 1);
    y -= std::lgammal(k + rho + 1);
    y += X.GetLogShape();
    return y;
}

template <typename IntType>
double YuleRand<IntType>::F(const IntType &k) const
{
    if (k < 1)
        return 0.0;
    double y = lgamma1pRo;
    y += math::lfact(k - 1);
    y -= std::lgammal(k + rho + 1);
    double logk = std::log(k);
    return -std::expm1(y + logk);
}

template <typename IntType>
double YuleRand<IntType>::S(const IntType &k) const
{
    if (k < 1)
        return 1.0;
    double y = lgamma1pRo;
    y += math::lfact(k - 1);
    y -= std::lgammal(k + rho + 1);
    y = std::exp(y);
    return k * y;
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType YuleRand<IntType>::operator()(G &g) const
{
    return variateGivenPareto(X(g), g);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType YuleRand<IntType>::Variate(double shape, G &g)
{
    if (!(shape > 0.0))
        throw std::invalid_argument("Yule distribution: shape should be positive");
    return variateGivenPareto(ParetoRand<double>::StandardVariate(shape, g), g);
}

template <typename IntType>
template <std::uniform_random_bit_generator G>
IntType YuleRand<IntType>::variateGivenPareto(double pareto, G &g)
{
    /// a success probability 1 / pareto that underflows leaves a geometric variate far beyond the integer range
    constexpr IntType largest = std::numeric_limits<IntType>::max();
    const double probability = 1.0 / pareto;
    if (!(probability > 0))
        return largest;
    const IntType failures = GeometricRand<IntType>::Variate(probability, g);
    return (failures < largest) ? failures + 1 : largest;
}

template <typename IntType>
long double YuleRand<IntType>::Mean() const
{
    return (rho <= 1) ? INFINITY : rho / (rho - 1);
}

template <typename IntType>
long double YuleRand<IntType>::Variance() const
{
    if (rho <= 2)
        return INFINITY;
    double aux = rho / (rho - 1);
    return aux * aux / (rho - 2);
}

template <typename IntType>
IntType YuleRand<IntType>::Mode() const
{
    return 1;
}

template <typename IntType>
long double YuleRand<IntType>::Skewness() const
{
    if (rho <= 3)
        return INFINITY;
    long double skewness = rho + 1;
    skewness *= skewness;
    skewness *= std::sqrt(rho - 2);
    return skewness / (rho * (rho - 3));
}

template <typename IntType>
long double YuleRand<IntType>::ExcessKurtosis() const
{
    if (rho <= 4)
        return INFINITY;
    long double numerator = 11 * rho * rho - 49;
    numerator *= rho;
    numerator -= 22;
    long double denominator = rho * (rho - 4) * (rho - 3);
    return rho + 3 + numerator / denominator;
}

} // namespace randlib

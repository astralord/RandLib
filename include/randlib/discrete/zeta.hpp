#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The ZetaRand class <BR>
 * Zeta distribution
 *
 * P(X = k) = 1 / (k^s * ζ(s))
 *
 * Notation: X ~ Zeta(s)
 */
template <typename IntType>
class ZetaRand final : public DiscreteDistribution<IntType>
{
    double s = 2; ///< exponent
    double sm1 = 1; ///< s - 1
    double zetaS = randlib::consts::pi_sq / 6.0; ///< ζ(s), where ζ stands for Riemann zeta-function
    double logZetaS = 2 * randlib::consts::lnpi - randlib::consts::ln2 - randlib::consts::ln3;///< ln(ζ(s))
    double b = 0.5; ///< 1 - 2^(1-s)

public:
    explicit ZetaRand(double exponent = 2.0);
    [[nodiscard]] std::string Name() const override;

    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] IntType MinValue() const override { return 1; }
    [[nodiscard]] IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetExponent(double exponent);
    [[nodiscard]] inline double GetExponent() const noexcept { return s; }

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
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
    [[nodiscard]] long double Moment(int n) const;
    [[nodiscard]] long double ThirdMoment() const override { return Moment(3); }
    [[nodiscard]] long double FourthMoment() const override { return Moment(4); }

    [[nodiscard]] inline long double GetZetaFunction() const noexcept { return zetaS; }
    [[nodiscard]] inline long double GetLogZetaFunction() const noexcept { return logZetaS; }
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>
#include <randlib/continuous/pareto.hpp>

namespace randlib
{

template < typename IntType >
ZetaRand<IntType>::ZetaRand(double exponent)
{
    SetExponent(exponent);
}

template < typename IntType >
std::string ZetaRand<IntType>::Name() const
{
    return "Zeta(" + this->toStringWithPrecision(GetExponent()) + ")";
}

template < typename IntType >
void ZetaRand<IntType>::SetExponent(double exponent)
{
    if (!(exponent > 1.0))
        throw std::invalid_argument("Zeta distribution: exponent should be greater than 1, but it's equal to " + std::to_string(exponent));
    s = exponent;
    sm1 = s - 1.0;
    zetaS = math::riemannZeta(s);
    logZetaS = std::log(zetaS);
    b = -std::expm1l(-sm1 * randlib::consts::ln2);
}

template < typename IntType >
double ZetaRand<IntType>::logP(const IntType & k) const
{
    return (k < 1) ? -INFINITY : -logZetaS - s * std::log(k);
}

template < typename IntType >
double ZetaRand<IntType>::F(const IntType & k) const
{
    return (k < 1) ? 0.0 : math::harmonicNumber(s, k) / zetaS;
}

template < typename IntType >
template <std::uniform_random_bit_generator G>
IntType ZetaRand<IntType>::operator()(G &g) const
{
    /// Luc Devroye, p. 551
    /// rejection sampling from rounded down Pareto distribution
    size_t iter = 0;
    do {
        float X = std::floor(ParetoRand<float>::StandardVariate(sm1, g));
        float V = detail::uniform01<float>(g);
        if (X < 1e4) {
            float T = std::pow(1.0 + 1.0 / X, sm1);
            /// there was a typo in the book - '<=' instead of '>'
            if (V * X * (T - 1) <= b * T)
                return X;
        }
        else {
            long double logT = sm1 * std::log1pl(1.0 / X);
            long double TM1 = std::expm1l(logT);
            if (V * X * TM1 <= b * std::exp(logT))
                return X < this->MaxValue() ? (IntType)X : this->MaxValue();
        }

    } while (++iter <= ProbabilityDistribution<IntType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Zeta distribution: sampling failed");
}

template < typename IntType >
long double ZetaRand<IntType>::Mean() const
{
    return (s > 2) ? math::riemannZeta(sm1) / zetaS : INFINITY;
}

template < typename IntType >
long double ZetaRand<IntType>::Variance() const
{
    if (s <= 3)
        return INFINITY;
    double y = Mean();
    double z = math::riemannZeta(s - 2) / zetaS;
    return z - y * y;
}

template < typename IntType >
IntType ZetaRand<IntType>::Mode() const
{
    return 1;
}

template < typename IntType >
long double ZetaRand<IntType>::Skewness() const
{
    if (s <= 4)
        return INFINITY;
    long double z1 = math::riemannZeta(sm1), z1Sq = z1 * z1;
    long double z2 = math::riemannZeta(s - 2);
    long double z3 = math::riemannZeta(s - 3);
    long double z = zetaS, zSq = z * z;
    long double logskew = zSq * z3;
    logskew -= 3 * z2 * z1 * z;
    logskew += 2 * z1 * z1Sq;
    logskew = std::log(logskew);
    logskew -= 1.5 * std::log(z * z2 - z1Sq);
    logskew -= 2 * logZetaS;
    return std::exp(logskew);
}

template < typename IntType >
long double ZetaRand<IntType>::ExcessKurtosis() const
{
    if (s <= 5)
        return INFINITY;
    long double mean = Mean();
    long double secondMoment = this->SecondMoment();
    long double thirdMoment = ThirdMoment();
    long double fourthMoment = FourthMoment();
    long double meanSq = mean * mean;
    long double variance = secondMoment - meanSq;
    long double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    long double denominator = variance * variance;
    return numerator / denominator - 3.0;
}

template < typename IntType >
long double ZetaRand<IntType>::Moment(int n) const
{
    return (s > n + 1) ? math::riemannZeta(s - n) / zetaS : INFINITY;
}

} // namespace randlib

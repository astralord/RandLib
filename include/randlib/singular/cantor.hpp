#pragma once

#include <array>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/singular.hpp>

namespace randlib
{

/**
 * @brief The CantorRand class <BR>
 * Cantor distribution
 *
 * Notation X ~ Cantor()
 */
class CantorRand final : public SingularDistribution
{
    static constexpr int n = 30;
    /// powers 3^-1, ..., 3^-n
    static constexpr std::array<double, n> table = [] {
        std::array<double, n> powers{};
        double power = 1.0;
        for (double &p : powers)
            p = power /= 3.0;
        return powers;
    }();

public:
    CantorRand();
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] double MinValue() const override { return 0; }
    [[nodiscard]] double MaxValue() const override { return 1; }

    [[nodiscard]] double F(const double & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] double operator()(G &g) const;
protected:
    double variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] double Median() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    double quantileImpl(double p, double initValue) const override;
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p, double initValue) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/discrete/bernoulli.hpp>
#include <randlib/continuous/uniform.hpp>

namespace randlib
{

inline CantorRand::CantorRand()
{
}

inline std::string CantorRand::Name() const
{
    return "Cantor";
}

inline double CantorRand::F(const double & x) const
{
    if (x <= 0.0)
        return 0.0;
    if (x >= 1.0)
        return 1.0;
    double a = 0.0, b = 1.0;
    double v = 0.5, d = 0.5, delta = 1.0;
    while (delta > MIN_POSITIVE)
    {
        delta = (b - a) / 3.0;
        if (x < a + delta)
        {
            b = a + delta;
            d *= 0.5;
            v -= d;
        }
        else if (x > b - delta)
        {
            a = b - delta;
            d *= 0.5;
            v += d;
        }
        else
            return v;
    }
    return v;
}

template <std::uniform_random_bit_generator G>
double CantorRand::operator()(G &g) const
{
    long double sum = 0.0;
    for (int i = 0; i != n; ++i) {
        sum += table[i] * BernoulliRand::StandardVariate(g);
    }
    return sum + sum;
}

inline long double CantorRand::Mean() const
{
    return 0.5;
}

inline long double CantorRand::Variance() const
{
    return 0.125;
}

inline double CantorRand::quantileImpl(double p, double) const
{
    return math::findRoot([this, p] (double x)
    {
        return F(x) - p;
    }, 0.0, 1.0);
}

inline double CantorRand::quantileImpl(double p) const
{
    return quantileImpl(p, p);
}

inline double CantorRand::quantileImpl1m(double p, double) const
{
    return math::findRoot([this, p] (double x)
    {
        return S(x) - p;
    }, 0.0, 1.0);
}

inline double CantorRand::quantileImpl1m(double p) const
{
    return quantileImpl1m(p, 1.0 - p);
}

inline std::complex<double> CantorRand::CFImpl(double t) const
{
    double prod = 1.0;
    for (int i = 0; i != n; ++i)
        prod *= std::cos(t * table[i]);
    std::complex<double> y(0.0, 0.5 * t);
    y = std::exp(y);
    return y * prod;
}

inline double CantorRand::Median() const
{
    return 1.0 / 3;
}

inline long double CantorRand::Skewness() const
{
    return 0.0l;
}

inline long double CantorRand::ExcessKurtosis() const
{
    return -1.6l;
}

} // namespace randlib

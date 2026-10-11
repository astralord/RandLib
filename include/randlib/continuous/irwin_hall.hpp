#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/uniform.hpp>

namespace randlib
{

/**
 * @brief The IrwinHallRand class <BR>
 * Irwin-Hall distribution
 *
 * f(x | n) = 0.5 / (n - 1)! * sum_{k=0}^n (-1)^k * C(n,k) * (x - k) ^ (n - 1) * sign(x - k)
 *
 * Notation: X ~ IH(n)
 *
 * Related distributions: <BR>
 * X ~ Y_1 + Y_2 + ... + Y_n, where Y_i ~ U(0,1)
 */
template <typename RealType>
class IrwinHallRand final : public ContinuousDistribution<RealType>
{
    int n = 1; ///< parameter of the distribution
    UniformRand<RealType> U{};

public:
    explicit IrwinHallRand(size_t number);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return n; }

    void SetNumber(int number);
    [[nodiscard]] inline int GetNumber() const noexcept { return n; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

namespace randlib
{

template < typename RealType >
IrwinHallRand<RealType>::IrwinHallRand(size_t number)
{
    SetNumber(number);
}

template < typename RealType >
std::string IrwinHallRand<RealType>::Name() const
{
    return "Irwin-Hall(" + this->toStringWithPrecision(GetNumber()) + ")";
}

template < typename RealType >
void IrwinHallRand<RealType>::SetNumber(int number)
{
    if (!(number > 0))
        throw std::invalid_argument("Irwin-Hall distribution: number should be positive");
    n = number;
}

template < typename RealType >
double IrwinHallRand<RealType>::f(const RealType & x) const
{
    if (x < 0 || x > n)
        return 0.0;

    /// Check simplest cases
    if (n == 1)
        return 1.0;
    if (n == 2)
        return (x <= 1.0) ? x : 2.0 - x;
    /// General case
    double sum = 0.0;
    int last = std::floor(x);
    for (int i = 0; i <= last; ++i) {
        double add = (n - 1) * std::log(x - i);
        add -= math::lfact(i);
        add -= math::lfact(n - i);
        add = std::exp(add);
        sum += (i % 2) ? -add : add;
    }
    return n * sum;
}

template < typename RealType >
double IrwinHallRand<RealType>::logf(const RealType & x) const
{
    return std::log(f(x));
}

template < typename RealType >
double IrwinHallRand<RealType>::F(const RealType & x) const
{
    if (x <= 0)
        return 0.0;
    if (x >= n)
        return 1.0;

    /// Check simplest cases
    if (n == 1)
        return x;
    if (n == 2) {
        if (x <= 1.0)
            return 0.5 * x * x;
        double temp = 2.0 - x;
        return 1.0 - 0.5 * temp * temp;
    }
    /// General case
    double sum = 0.0;
    int last = std::floor(x);
    for (int i = 0; i <= last; ++i) {
        double add = n * std::log(x - i);
        add -= math::lfact(i);
        add -= math::lfact(n - i);
        add = std::exp(add);
        sum += (i % 2) ? -add : add;
    }
    return sum;
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType IrwinHallRand<RealType>::operator()(G &g) const
{
    RealType sum = 0.0;
    for (int i = 0; i != n; ++i)
        sum += U(g);
    return sum;
}



template < typename RealType >
long double IrwinHallRand<RealType>::Mean() const
{
    return 0.5 * n;
}

template < typename RealType >
long double IrwinHallRand<RealType>::Variance() const
{
    return n / 12.0l;
}

template < typename RealType >
std::complex<double> IrwinHallRand<RealType>::CFImpl(double t) const
{
    return std::pow(U.CF(t), n);
}

template < typename RealType >
RealType IrwinHallRand<RealType>::Median() const
{
    return 0.5 * n;
}

template < typename RealType >
RealType IrwinHallRand<RealType>::Mode() const
{
    return 0.5 * n;
}

template < typename RealType >
long double IrwinHallRand<RealType>::Skewness() const
{
    return 0.0l;
}

template < typename RealType >
long double IrwinHallRand<RealType>::ExcessKurtosis() const
{
    return -1.2 / n;
}

} // namespace randlib

#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The TriangularRand class <BR>
 * Triangular distribution
 *
 * Notation: X ~ Tri(a, b, c)
 */
template <typename RealType>
class TriangularRand final : public ContinuousDistribution<RealType>
{
    double a = 0; ///< min value
    double b = 2; ///< max value
    double c = 1; ///< mode
    double constForGenerator = 1; ///< (c - a) / (b - a)
    double coefGenerator1 = 1; ///< (b - a) * (c - a)
    double coefGenerator2 = 1; ///< (b - a) * (b - c)

    void SetConstantsForGenerator();

public:
    TriangularRand(double lowerLimit = 0, double mode = 0.5, double upperLimit = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return a; }
    [[nodiscard]] RealType MaxValue() const override { return b; }

    void SetParameters(double lowerLimit, double mode, double upperLimit);

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:

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

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template < typename RealType >
TriangularRand<RealType>::TriangularRand(double lowerLimit, double mode, double upperLimit)
{
    SetParameters(lowerLimit, mode, upperLimit);
}

template < typename RealType >
std::string TriangularRand<RealType>::Name() const
{
    return "Triangular("
            + this->toStringWithPrecision(MinValue()) + ", "
            + this->toStringWithPrecision(Mode()) + ", "
            + this->toStringWithPrecision(MaxValue()) + ")";
}

template < typename RealType >
void TriangularRand<RealType>::SetParameters(double lowerLimit, double mode, double upperLimit)
{
    if (!(lowerLimit < upperLimit))
        throw std::invalid_argument("Triangular distribution: lower limit should be smaller than upper limit");
    if (!(lowerLimit <= mode && mode <= upperLimit))
        throw std::invalid_argument("Triangular distribution: mode should lie between the lower and upper limits");
    a = lowerLimit;
    c = mode;
    b = upperLimit;
    SetConstantsForGenerator();
}

template < typename RealType >
void TriangularRand<RealType>::SetConstantsForGenerator()
{
    constForGenerator = (c - a) / (b - a);
    coefGenerator1 = (b - a) * (c - a);
    coefGenerator2 = (b - a) * (b - c);
}

template < typename RealType >
double TriangularRand<RealType>::f(const RealType & x) const
{
    if (x < a || x > b)
        return 0;
    if (x < c)
        return 2.0 * (x - a) / coefGenerator1;
    if (x == c)
        return 2.0 / (b - a);
    if (x < b)
        return 2.0 * (b - x) / coefGenerator2;
    return 0;
}

template < typename RealType >
double TriangularRand<RealType>::logf(const RealType & x) const
{
    return std::log(f(x));
}

template < typename RealType >
double TriangularRand<RealType>::F(const RealType & x) const
{
    if (x <= a)
        return 0.0;
    if (x <= c)
        return (x - a) * (x - a) / coefGenerator1;
    if (x < b)
        return 1.0 - (b - x) * (b - x) / coefGenerator2;
    return 1.0;
}

template < typename RealType >
double TriangularRand<RealType>::S(const RealType & x) const
{
    if (x <= a)
        return 1.0;
    if (x <= c)
        return 1.0 - (x - a) * (x - a) / coefGenerator1;
    if (x < b)
        return (b - x) * (b - x) / coefGenerator2;
    return 0.0;
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType TriangularRand<RealType>::operator()(G &g) const
{
    RealType U = detail::uniform01<RealType>(g);
    if (U < constForGenerator)
        return a + std::sqrt(U * coefGenerator1);
    return b - std::sqrt((1 - U) * coefGenerator2);
}

template < typename RealType >
long double TriangularRand<RealType>::Mean() const
{
    return (a + b + c) / 3.0;
}

template < typename RealType >
long double TriangularRand<RealType>::Variance() const
{
    return (a * (a - b) + b * (b - c) + c * (c - a)) / 18.0;
}

template < typename RealType >
std::complex<double> TriangularRand<RealType>::CFImpl(double t) const
{
    if (c == a || c == b) {
        /// Right triangle: for c = a, φ(t) = e^{ita} 2/h² ((1 - e^{ith}) / t² + ih / t) with h = b - a;
        /// for c = b, reflect X = a + b - Y
        double h = b - a, s = (c == a) ? t : -t;
        std::complex<double> eish(std::cos(s * h), std::sin(s * h));
        std::complex<double> y = 2.0 / (h * h) * ((1.0 - eish) / (s * s) + std::complex<double>(0.0, h / s));
        double shift = (c == a) ? a * t : b * t;
        return std::complex<double>(std::cos(shift), std::sin(shift)) * y;
    }
    double bmc = b - c, bma = b - a, cma = c - a;
    double at = a * t, bt = b * t, ct = c * t;
    std::complex<double> x(bmc * std::cos(at), bmc * std::sin(at));
    std::complex<double> y(bma * std::cos(ct), bma * std::sin(ct));
    std::complex<double> z(cma * std::cos(bt), cma * std::sin(bt));
    std::complex<double> numerator = x - y + z;
    /// in order to avoid numerical errors
    if (t < 1e-10 && std::fabs(numerator.real()) < 1e-10)
        return 1;
    double denominator = bma * cma * bmc * t * t;
    std::complex<double> frac = -numerator / denominator;
    return frac + frac;
}

template < typename RealType >
RealType TriangularRand<RealType>::Median() const
{
    if (c + c > a + b)
        return a + std::sqrt(0.5 * coefGenerator1);
    return b - std::sqrt(0.5 * coefGenerator2);
}

template < typename RealType >
RealType TriangularRand<RealType>::Mode() const
{
    return c;
}

template < typename RealType >
long double TriangularRand<RealType>::Skewness() const
{
    double numerator = randlib::consts::sqrt2;
    numerator *= (a + b - c - c);
    numerator *= (a + a - b - c);
    numerator *= (a - b - b + c);
    double denominator = a * (a - b);
    denominator += b * (b - c);
    denominator += c * (c - a);
    denominator *= std::sqrt(denominator);
    return 0.2 * numerator / denominator;
}

template < typename RealType >
long double TriangularRand<RealType>::ExcessKurtosis() const
{
    return -0.6;
}

} // namespace randlib

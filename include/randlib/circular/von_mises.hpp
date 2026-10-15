#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/circular.hpp>

namespace randlib
{

/**
 * @brief The VonMisesRand class <BR>
 * Von-Mises distribution
 *
 * Notation: X ~ Von-Mises(μ, k)
 */
template <typename RealType>
class VonMisesRand final : public CircularDistribution<RealType>
{
    double k = 1; ///< concentration
    double logI0k = 1; ///< log(I_0(k)), where I_0 stands for modified Bessel function of the first kind of order 0
    double s = randlib::consts::pi / randlib::consts::e; ///< generator coefficient
    int p = 12; /// coefficient for faster cdf calculation
    static constexpr double CK = 10.5;

public:
    VonMisesRand(double location = 0, double concentration = 1);

    [[nodiscard]] std::string Name() const override;

    void SetConcentration(double concentration);
    [[nodiscard]] inline double GetConcentration() const noexcept { return k; }

private:
    double cdfSeries(double x) const;
    double cdfErfcAux(RealType x) const;
    double cdfErfc(RealType x) const;
    double ccdfErfc(RealType x) const;

public:
    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

public:
    [[nodiscard]] long double CircularMean() const override;
    [[nodiscard]] long double CircularVariance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template <typename RealType>
VonMisesRand<RealType>::VonMisesRand(double location, double concentration) : CircularDistribution<RealType>(location)
{
    SetConcentration(concentration);
}

template <typename RealType>
std::string VonMisesRand<RealType>::Name() const
{
    return "von Mises(" + this->toStringWithPrecision(this->GetLocation()) + ", " +
           this->toStringWithPrecision(this->GetConcentration()) + ")";
}

template <typename RealType>
void VonMisesRand<RealType>::SetConcentration(double concentration)
{
    if (!(concentration > 0.0))
        throw std::invalid_argument("von Mises distribution: concentration parameter should be positive");
    k = concentration;
    logI0k = math::logBesselI(0, k);
    s = (k > 1.3) ? 1.0 / std::sqrt(k) : randlib::consts::pi * std::exp(-k);
    if (k < CK) {
        /// set coefficients for cdf
        static constexpr double a[] = {28.0, 0.5, 100.0, 5.0};
        p = std::ceil(a[0] + a[1] * k - a[2] / (k + a[3]));
    }
}

template <typename RealType>
double VonMisesRand<RealType>::cdfSeries(double x) const
{
    /// backwards recursion
    double sinX = std::sin(x), cosX = std::cos(x);
    double px = p * x;
    double sn = std::sin(px), cn = std::cos(px);
    double R = 0, V = 0;
    for (int n = p - 1; n > 0; --n) {
        double temp = sn;
        sn *= cosX;
        sn -= cn * sinX;
        cn *= cosX;
        cn += temp * sinX;
        R += 2 * n / k;
        R = 1.0 / R;
        V += sn / n;
        V *= R;
    }
    V /= randlib::consts::pi;
    V += 0.5 * x / randlib::consts::pi;
    V += 0.5;
    return V;
}

template <typename RealType>
double VonMisesRand<RealType>::cdfErfcAux(RealType x) const
{
    /// Normal cdf approximation
    double c = 24.0 * k;
    double v = c - 56.0;
    double r = std::sqrt((54.0 / (347.0 / v + 26.0 - c) - 6.0 + c) / 12.0);
    double z = -std::sin(0.5 * x) * r;
    double twoZSq = 2.0 * z * z;
    v -= twoZSq - 3.0;
    double y = (c - 2 * twoZSq - 16.0) / 3.0;
    y = ((twoZSq + 1.75) * twoZSq + 83.5) / v - y;
    y *= y;
    return z * (1.0 - twoZSq / y);
}

template <typename RealType>
double VonMisesRand<RealType>::cdfErfc(RealType x) const
{
    double arg = cdfErfcAux(x);
    return 0.5 * std::erfc(arg);
}

template <typename RealType>
double VonMisesRand<RealType>::ccdfErfc(RealType x) const
{
    double arg = cdfErfcAux(x);
    return 0.5 * std::erfc(-arg);
}

template <typename RealType>
double VonMisesRand<RealType>::f(const RealType &x) const
{
    return (x < this->MinValue() || x > this->MaxValue()) ? 0.0 : std::exp(logf(x));
}

template <typename RealType>
double VonMisesRand<RealType>::logf(const RealType &x) const
{
    if (x < this->MinValue() || x > this->MaxValue())
        return -INFINITY;
    double y = k * std::cos(x - this->loc) - logI0k;
    y -= randlib::consts::lnpi + randlib::consts::ln2;
    return y;
}

template <typename RealType>
double VonMisesRand<RealType>::F(const RealType &x) const
{
    if (x <= this->MinValue())
        return 0.0;
    if (x >= this->MaxValue())
        return 1.0;
    double xAdj = x - this->loc;
    return (k < CK) ? cdfSeries(xAdj) : cdfErfc(xAdj);
}

template <typename RealType>
double VonMisesRand<RealType>::S(const RealType &x) const
{
    if (x <= this->MinValue())
        return 1.0;
    if (x >= this->MaxValue())
        return 0.0;
    double xAdj = x - this->loc;
    return (k < CK) ? 1.0 - cdfSeries(xAdj) : ccdfErfc(xAdj);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType VonMisesRand<RealType>::operator()(G &g) const
{
    /// Generating von Mises variates by the ratio-of-uniforms method
    /// Lucio Barabesi. Dipartimento di Metodi Quantitativi, Universiteta di Siena
    size_t iter = 0;
    do {
        RealType U = detail::uniform01<RealType>(g);
        RealType V = 2 * detail::uniform01<RealType>(g) - 1;
        RealType theta = s * V / U;
        if ((std::fabs(theta) <= randlib::consts::pi) &&
            ((k * theta * theta < 4.0 * (1.0 - U)) || (k * std::cos(theta) >= 2 * std::log(U) + k)))
            return this->loc + theta;
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("von Mises distribution: sampling failed");
}

template <typename RealType>
long double VonMisesRand<RealType>::CircularMean() const
{
    return this->loc;
}

template <typename RealType>
long double VonMisesRand<RealType>::CircularVariance() const
{
    double var = math::logBesselI(1, k);
    var -= logI0k;
    return -std::expm1l(var);
}

template <typename RealType>
std::complex<double> VonMisesRand<RealType>::CFImpl(double t) const
{
    /// I_t(κ) / I_0(κ) is the trigonometric moment, which equals E[exp(itX)] only for integer t
    if (t != std::round(t))
        return CircularDistribution<RealType>::CFImpl(t);
    double tloc = t * this->loc;
    double cosTloc = std::cos(tloc), sinTloc = std::sin(tloc);
    std::complex<double> y(cosTloc, sinTloc);
    double z = math::logBesselI(t, k) - logI0k;
    return y * std::exp(z);
}

template <typename RealType>
RealType VonMisesRand<RealType>::Median() const
{
    return this->loc;
}

template <typename RealType>
RealType VonMisesRand<RealType>::Mode() const
{
    return this->loc;
}

} // namespace randlib

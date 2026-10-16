#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/discrete/poisson.hpp>

namespace randlib
{

/**
 * @brief The NoncentralChiSquaredRand class <BR>
 * Noncentral Chi-Squared distribution
 *
 * Notation: X ~ χ'^2(k, λ)
 *
 * Related distributions: <BR>
 * If X ~ χ'^2(k, 0), then X ~ χ^2(k) <BR>
 * X ~ χ^2(k + 2J), where J ~ Po(λ)
 */
template <typename RealType>
class NoncentralChiSquaredRand final : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    double k = 1; ///< degree
    double lambda = 2; ///< noncentrality λ
    double halfK = 0.5; ///< k / 2
    double halfLambda = 1; ///< λ / 2
    double sqrtLambda = randlib::consts::sqrt2; ///< √λ
    double logLambda = randlib::consts::ln2; ///< log(λ)

    PoissonRand<int> Y{};

public:
    explicit NoncentralChiSquaredRand(double degree = 1, double noncentrality = 0);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double degree, double noncentrality);
    [[nodiscard]] inline double GetDegree() const noexcept { return k; }
    [[nodiscard]] inline double GetNoncentrality() const noexcept { return lambda; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;

private:
    template <std::uniform_random_bit_generator G>
    RealType variateForDegreeEqualOne(G &g) const;

public:
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType Variate(double degree, double noncentrality, G &g);
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/gamma.hpp>
#include <randlib/continuous/normal.hpp>

namespace randlib
{

template <typename RealType>
NoncentralChiSquaredRand<RealType>::NoncentralChiSquaredRand(double degree, double noncentrality)
{
    SetParameters(degree, noncentrality);
}

template <typename RealType>
std::string NoncentralChiSquaredRand<RealType>::Name() const
{
    return "Noncentral Chi-Squared(" + this->toStringWithPrecision(GetDegree()) + ", " +
           this->toStringWithPrecision(GetNoncentrality()) + ")";
}

template <typename RealType>
void NoncentralChiSquaredRand<RealType>::SetParameters(double degree, double noncentrality)
{
    if (!(degree > 0.0))
        throw std::invalid_argument("Noncentral Chi-Squared distribution: degree parameter should be positive");
    if (!(noncentrality > 0.0))
        throw std::invalid_argument("Noncentral Chi-Squared distribution: noncentrality parameter should be positive");

    k = degree;
    halfK = 0.5 * k;

    lambda = noncentrality;
    halfLambda = 0.5 * lambda;
    sqrtLambda = std::sqrt(lambda);
    logLambda = std::log(lambda);

    if (k < 1)
        Y.SetRate(halfLambda);
}

template <typename RealType>
double NoncentralChiSquaredRand<RealType>::f(const RealType &x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (k == 2)
            return 0.5 * std::exp(-halfLambda);
        return (k > 2) ? 0.0 : INFINITY;
    }
    return std::exp(logf(x));
}

template <typename RealType>
double NoncentralChiSquaredRand<RealType>::logf(const RealType &x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0) {
        if (k == 2)
            return -halfLambda - randlib::consts::ln2;
        return (k > 2) ? -INFINITY : INFINITY;
    }
    double halfKm1 = halfK - 1;
    double y = math::logBesselI(halfKm1, std::sqrt(lambda * x));
    double z = halfKm1 * (std::log(x) - logLambda);
    z -= x + lambda;
    return y + 0.5 * z - randlib::consts::ln2;
}

template <typename RealType>
double NoncentralChiSquaredRand<RealType>::F(const RealType &x) const
{
    if (x <= 0.0)
        return 0.0;
    if (std::isinf(x))
        return 1.0;
    double halfX = 0.5 * x;
    double sqrtHalfX = std::sqrt(halfX);
    double logHalfX = std::log(halfX);
    return math::MarcumP(halfK, halfLambda, halfX, sqrtLambda / randlib::consts::sqrt2, sqrtHalfX,
                         logLambda - randlib::consts::ln2, logHalfX);
}

template <typename RealType>
double NoncentralChiSquaredRand<RealType>::S(const RealType &x) const
{
    if (x <= 0.0)
        return 1.0;
    if (std::isinf(x))
        return 0.0;
    double halfX = 0.5 * x;
    double sqrtHalfX = std::sqrt(halfX);
    double logHalfX = std::log(halfX);
    return math::MarcumQ(halfK, halfLambda, halfX, sqrtLambda / randlib::consts::sqrt2, sqrtHalfX,
                         logLambda - randlib::consts::ln2, logHalfX);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType NoncentralChiSquaredRand<RealType>::variateForDegreeEqualOne(G &g) const
{
    RealType y = sqrtLambda + detail::standardNormal<RealType>(g);
    return y * y;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType NoncentralChiSquaredRand<RealType>::Variate(double degree, double noncentrality, G &g)
{
    if (!(degree > 0.0))
        throw std::invalid_argument("Noncentral Chi-Squared distribution: degree parameter should be positive");
    if (!(noncentrality > 0.0))
        throw std::invalid_argument("Noncentral Chi-Squared distribution: noncentrality parameter should be positive");

    if (degree >= 1) {
        RealType rv = (degree == 1) ? 0.0 : 2 * GammaDistribution<RealType>::StandardVariate(0.5 * degree - 0.5, g);
        RealType y = std::sqrt(noncentrality) + detail::standardNormal<RealType>(g);
        return rv + y * y;
    }
    RealType shape = 0.5 * degree + PoissonRand<int>::Variate(0.5 * noncentrality, g);
    return 2 * GammaDistribution<RealType>::StandardVariate(shape, g);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType NoncentralChiSquaredRand<RealType>::operator()(G &g) const
{
    if (k < 1)
        return 2 * GammaDistribution<RealType>::StandardVariate(halfK + Y(g), g);
    double X = variateForDegreeEqualOne(g);
    if (k > 1)
        X += 2 * GammaDistribution<RealType>::StandardVariate(halfK - 0.5, g);
    return X;
}

template <typename RealType>
void NoncentralChiSquaredRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    if (k >= 1) {
        for (RealType &var : outputData)
            var = variateForDegreeEqualOne(g);
        double halfKmHalf = halfK - 0.5;
        if (halfKmHalf == 0)
            return;
        for (RealType &var : outputData)
            var += 2 * GammaDistribution<RealType>::StandardVariate(halfKmHalf, g);
    } else {
        for (RealType &var : outputData)
            var = 2 * GammaDistribution<RealType>::StandardVariate(halfK + Y(g), g);
    }
}

template <typename RealType>
long double NoncentralChiSquaredRand<RealType>::Mean() const
{
    return k + lambda;
}

template <typename RealType>
long double NoncentralChiSquaredRand<RealType>::Variance() const
{
    return 2 * (k + 2 * lambda);
}

template <typename RealType>
RealType NoncentralChiSquaredRand<RealType>::Mode() const
{
    /// For k < 2 the density is unbounded at 0; for k = 2 it is finite there, and the
    /// maximum moves away from 0 for large enough noncentrality
    if (k < 2)
        return 0.0;
    RealType mode = ContinuousDistribution<RealType>::Mode();
    return (mode > 0 && this->f(mode) > this->f(0)) ? mode : 0.0;
}

template <typename RealType>
long double NoncentralChiSquaredRand<RealType>::Skewness() const
{
    long double y = k + 2 * lambda;
    y = 2.0 / y;
    long double z = y * std::sqrt(y);
    return z * (k + 3 * lambda);
}

template <typename RealType>
long double NoncentralChiSquaredRand<RealType>::ExcessKurtosis() const
{
    long double y = k + 2 * lambda;
    return 12 * (k + 4 * lambda) / (y * y);
}

template <typename RealType>
std::complex<double> NoncentralChiSquaredRand<RealType>::CFImpl(double t) const
{
    std::complex<double> aux(1, -2 * t);
    std::complex<double> y(0, lambda * t);
    y /= aux;
    y -= halfK * std::log(aux);
    return std::exp(y);
}

} // namespace randlib

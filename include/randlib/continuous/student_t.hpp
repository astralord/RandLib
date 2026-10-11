#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/nakagami.hpp>

namespace randlib
{

/**
 * @brief The StudentTRand class <BR>
 * Student's t-distribution
 *
 * Notation: X ~ t(ν, μ, σ)
 * If X ~ t(1, μ, σ), then X ~ Cauchy(μ, σ)
 * X -> Normal(μ, σ) for t -> ∞
 */
template <typename RealType>
class StudentTRand final : public ContinuousDistribution<RealType>
{
    double nu = 1; ///< degree ν
    double mu = 0; ///< location μ
    double sigma = 1; ///< scale σ
    double logSigma = 0; ///< log(σ)
    NakagamiRand<RealType> Y{};
    double pdfCoef = -randlib::consts::lnpi; ///< coefficient for faster pdf calculation
    double nup1Half = 1; ///< 0.5 * (ν + 1)
    double logBetaFun = randlib::consts::lnpi; ///< log(B(0.5 * ν, 0.5))

public:
    explicit StudentTRand(double degree = 1.0, double location = 0.0, double scale = 1.0);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] RealType MinValue() const override { return -INFINITY; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetDegree(double degree);
    void SetLocation(double location);
    void SetScale(double scale);
    [[nodiscard]] inline double GetDegree() const noexcept { return nu; }
    [[nodiscard]] inline double GetLocation() const noexcept { return mu; }
    [[nodiscard]] inline double GetScale() const noexcept { return sigma; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;

private:
    double log1pSquareOverDegree(double x0) const;
    double lowerTail(double x0) const;

public:
    [[nodiscard]] double F(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:
    protected:
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/cauchy.hpp>

namespace randlib
{

template < typename RealType >
StudentTRand<RealType>::StudentTRand(double degree, double location, double scale)
{
    SetDegree(degree);
    SetLocation(location);
    SetScale(scale);
}

template < typename RealType >
std::string StudentTRand<RealType>::Name() const
{
    if (mu == 0.0 && sigma == 1.0)
        return "Student-t(" + this->toStringWithPrecision(GetDegree()) + ")";
    return "Student-t(" + this->toStringWithPrecision(GetDegree()) + ", "
                          + this->toStringWithPrecision(GetLocation()) + ", "
                          + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void StudentTRand<RealType>::SetDegree(double degree)
{
    if (!(degree > 0.0))
        throw std::invalid_argument("Student-t distribution: degree parameter should be positive");
    nu = degree;
    Y.SetParameters(0.5 * nu, 1.0);

    nup1Half = 0.5 * (nu + 1);
    pdfCoef = Y.GetLogGammaShapeRatio();
    pdfCoef -= 0.5 * randlib::consts::lnpi;
    logBetaFun = -pdfCoef;
    pdfCoef -= 0.5 * std::log(nu);
}

template < typename RealType >
void StudentTRand<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Student-t distribution: location should be finite");
    mu = location;
}

template < typename RealType >
void StudentTRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Student-t distribution: scale should be positive");
    sigma = scale;
    logSigma = std::log(sigma);
}

template < typename RealType >
double StudentTRand<RealType>::f(const RealType & x) const
{
    /// adjustment
    double x0 = x - mu;
    x0 /= sigma;
    double xSq = x0 * x0;
    if (nu == 1) /// Cauchy distribution
        return randlib::consts::inv_pi / (sigma * (1 + xSq));
    if (nu == 2) {
        double tmp = 2.0 + xSq;
        return 1.0 / (sigma * tmp * std::sqrt(tmp));
    }
    if (nu == 3) {
        double y = 3 + xSq;
        return 6 * randlib::consts::sqrt3 * randlib::consts::inv_pi / (sigma * y * y);
    }
    return std::exp(logf(x));
}

template < typename RealType >
double StudentTRand<RealType>::log1pSquareOverDegree(double x0) const
{
    /// log(1 + x0² / ν), also when x0² overflows
    double absX = std::fabs(x0);
    if (absX < 1e100)
        return std::log1pl(x0 * x0 / nu);
    double logRatio = 2 * std::log(absX) - std::log(nu);
    return logRatio + std::log1pl(std::exp(-logRatio));
}

template < typename RealType >
double StudentTRand<RealType>::logf(const RealType & x) const
{
    /// adjustment
    double x0 = x - mu;
    x0 /= sigma;
    double log1pXsqNu = log1pSquareOverDegree(x0);
    if (nu == 1) /// Cauchy distribution
        return -logSigma - randlib::consts::lnpi - log1pXsqNu;
    if (nu == 2)
        return -logSigma - 1.5 * (randlib::consts::ln2 + log1pXsqNu);
    if (nu == 3) {
        double y = -2 * (randlib::consts::ln3 + log1pXsqNu);
        y += randlib::consts::ln2 + 1.5 * randlib::consts::ln3 - randlib::consts::lnpi - logSigma;
        return y;
    }
    double y = -nup1Half * log1pXsqNu;
    return pdfCoef + y - logSigma;
}

template < typename RealType >
double StudentTRand<RealType>::lowerTail(double x0) const
{
    /// P(T < x0) for x0 < 0: I_t(ν/2, 1/2) / 2 with t = ν / (ν + x0²)
    if (nu == 1)
        return randlib::consts::inv_pi * ((x0 < -1.0) ? std::atan(-1.0 / x0) : randlib::consts::half_pi + std::atan(x0));
    if (x0 * x0 < nu) {
        /// Near the median t = ν / (ν + x0²) rounds to 1; use 1/2 - I_{x0²/(ν+x0²)}(1/2, ν/2) / 2 instead
        double xSq = x0 * x0, u = xSq / (nu + xSq);
        return 0.5 - 0.5 * math::ibeta(u, 0.5, 0.5 * nu, logBetaFun, std::log(xSq) - std::log(nu + xSq), std::log(nu) - std::log(nu + xSq));
    }
    double logT = -log1pSquareOverDegree(x0);
    double t = std::exp(logT);
    if (t < 1e-280) {
        /// t underflows: leading term t^(ν/2) / ((ν/2) B(ν/2, 1/2)) of the incomplete beta function
        return 0.5 * std::exp(0.5 * nu * logT - std::log(0.5 * nu) - logBetaFun);
    }
    double log1mT = std::log1pl(-t);
    return 0.5 * math::ibeta(t, 0.5 * nu, 0.5, logBetaFun, logT, log1mT);
}

template < typename RealType >
double StudentTRand<RealType>::F(const RealType & x) const
{
    double x0 = (x - mu) / sigma;
    if (x0 == 0.0)
        return 0.5;
    return (x0 < 0.0) ? lowerTail(x0) : 1.0 - lowerTail(-x0);
}

template < typename RealType >
double StudentTRand<RealType>::S(const RealType & x) const
{
    double x0 = (x - mu) / sigma;
    if (x0 == 0.0)
        return 0.5;
    return (x0 > 0.0) ? lowerTail(-x0) : 1.0 - lowerTail(x0);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType StudentTRand<RealType>::operator()(G &g) const
{
    if (nu == 1)
        return mu + sigma * CauchyRand<RealType>::StandardVariate(g);
    return mu + sigma * detail::standardNormal<RealType>(g) / Y(g);
}

template < typename RealType >
void StudentTRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    if (nu == 1) {
        for (RealType &var : outputData)
            var = mu + sigma * CauchyRand<RealType>::StandardVariate(g);
    }
    else {
        Y.Sample(g, outputData);
        for (RealType &var : outputData)
            var = mu + sigma * detail::standardNormal<RealType>(g) / var;
    }
}



template < typename RealType >
long double StudentTRand<RealType>::Mean() const
{
    return (nu > 1) ? mu : NAN;
}

template < typename RealType >
long double StudentTRand<RealType>::Variance() const
{
    if (nu > 2)
        return sigma * sigma * nu / (nu - 2);
    return (nu > 1) ? INFINITY : NAN;
}

template < typename RealType >
std::complex<double> StudentTRand<RealType>::CFImpl(double t) const
{
    double x = std::sqrt(nu) * t * sigma; // value of sqrt(nu) can be hashed
    double vHalf = 0.5 * nu;
    double y = vHalf * std::log(x);
    y -= Y.GetLogGammaFunction();
    y -= (vHalf - 1) * randlib::consts::ln2;
    y += math::logBesselK(vHalf, x);
    double costmu = std::cos(t * mu), sintmu = std::sin(t * mu);
    std::complex<double> cf(costmu, sintmu);
    return std::exp(y) * cf;
}

template < typename RealType >
RealType StudentTRand<RealType>::quantileImpl(double p) const
{
    double temp = p - 0.5;
    if (nu == 1)
        return std::tan(randlib::consts::pi * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (nu == 2)
        return sigma * 2.0 * temp * std::sqrt(0.5 / pq) + mu;
    if (nu == 4)
    {
        double alpha = 2 * std::sqrt(pq);
        double beta = std::cos(std::acos(alpha) / 3.0) / alpha - 1;
        return mu + sigma * 2 * math::sign(temp) * std::sqrt(beta);
    }
    return ContinuousDistribution<RealType>::quantileImpl(p);
}

template < typename RealType >
RealType StudentTRand<RealType>::quantileImpl1m(double p) const
{
    double temp = 0.5 - p;
    if (nu == 1)
        return std::tan(randlib::consts::pi * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (nu == 2)
        return sigma * 2 * temp * std::sqrt(0.5 / pq) + mu;
    if (nu == 4)
    {
        double alpha = 2 * std::sqrt(pq);
        double beta = std::cos(std::acos(alpha) / 3.0) / alpha - 1;
        return mu + sigma * 2 * math::sign(temp) * std::sqrt(beta);
    }
    return ContinuousDistribution<RealType>::quantileImpl1m(p);
}

template < typename RealType >
RealType StudentTRand<RealType>::Median() const
{
    return mu;
}

template < typename RealType >
RealType StudentTRand<RealType>::Mode() const
{
    return mu;
}

template < typename RealType >
long double StudentTRand<RealType>::Skewness() const
{
    return (nu > 3) ? 0.0 : NAN;
}

template < typename RealType >
long double StudentTRand<RealType>::ExcessKurtosis() const
{
    if (nu > 4)
        return 6.0 / (nu - 4);
    return (nu > 2) ? INFINITY : NAN;
}

} // namespace randlib

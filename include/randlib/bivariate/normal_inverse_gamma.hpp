#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous_bivariate.hpp>
#include <randlib/continuous/inverse_gamma.hpp>
#include <randlib/continuous/student_t.hpp>

#include <algorithm>

namespace randlib
{

/**
 * @brief The NormalInverseGammaRand class <BR>
 * Normal-inverse-gamma distribution
 *
 * Notation: X ~ NIG(μ, λ, α, β)
 */
template <typename RealType>
class NormalInverseGammaRand final
    : public ContinuousBivariateDistribution<StudentTRand<RealType>, InverseGammaRand<RealType>, RealType>
{
    Pair<RealType> variate(detail::EngineRef g) const override { return (*this)(g); }

    double mu = 0; ///< location μ
    double lambda = 1; ///< precision λ
    double alpha = 1; ///< shape α
    double beta = 1; ///< rate β
    double pdfCoef = 0.5 * (randlib::consts::lnpi - randlib::consts::ln2); ///< coefficient for faster pdf calculation

public:
    NormalInverseGammaRand(double location = 0, double precision = 1, double shape = 1, double rate = 1);
    [[nodiscard]] std::string Name() const override;

    void SetParameters(double location, double precision, double shape, double rate);
    [[nodiscard]] inline double GetLocation() const noexcept { return mu; }
    [[nodiscard]] inline double GetPrecision() const noexcept { return lambda; }
    [[nodiscard]] inline double GetShape() const noexcept { return alpha; }
    [[nodiscard]] inline double GetRate() const noexcept { return beta; }

    [[nodiscard]] double f(const Pair<RealType> &point) const override;
    [[nodiscard]] double logf(const Pair<RealType> &point) const override;
    [[nodiscard]] double F(const Pair<RealType> &point) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] Pair<RealType> operator()(G &g) const;

    long double Correlation() const override;

    [[nodiscard]] Pair<RealType> Mode() const override;
};

} // namespace randlib

#include <randlib/continuous/normal.hpp>

namespace randlib
{

template <typename RealType>
NormalInverseGammaRand<RealType>::NormalInverseGammaRand(double location, double precision, double shape, double rate)
{
    SetParameters(location, precision, shape, rate);
}

template <typename RealType>
std::string NormalInverseGammaRand<RealType>::Name() const
{
    return "Normal-Inverse-Gamma(" + this->toStringWithPrecision(GetLocation()) + ", " +
           this->toStringWithPrecision(GetPrecision()) + ", " + this->toStringWithPrecision(GetShape()) + ", " +
           this->toStringWithPrecision(GetRate()) + ")";
}

template <typename RealType>
void NormalInverseGammaRand<RealType>::SetParameters(double location, double precision, double shape, double rate)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Normal-Inverse-Gamma distribution: location should be finite");
    if (!(precision > 0.0))
        throw std::invalid_argument(
            "Precision of Normal-Inverse-Gamma distribution should be positive, but it's equal to " +
            std::to_string(precision));
    if (!(shape > 0.0))
        throw std::invalid_argument(
            "Shape of Normal-Inverse-Gamma distribution should be positive, but it's equal to " +
            std::to_string(shape));
    if (!(rate > 0.0))
        throw std::invalid_argument("Rate of Normal-Inverse-Gamma distribution should be positive, but it's equal to " +
                                    std::to_string(rate));

    mu = location;
    lambda = precision;

    this->Y.SetParameters(shape, rate);
    alpha = this->Y.GetShape();
    beta = this->Y.GetRate();
    this->X.SetDegree(2 * alpha);
    this->X.SetLocation(mu);
    this->X.SetScale(std::sqrt(beta / (alpha * lambda)));

    pdfCoef = 0.5 * std::log(0.5 * lambda / randlib::consts::pi);
    pdfCoef += alpha * this->Y.GetLogRate() - this->Y.GetLogGammaShape();
}

template <typename RealType>
double NormalInverseGammaRand<RealType>::f(const Pair<RealType> &point) const
{
    return (point.second > 0.0) ? std::exp(logf(point)) : 0.0;
}

template <typename RealType>
double NormalInverseGammaRand<RealType>::logf(const Pair<RealType> &point) const
{
    RealType sigmaSq = point.second;
    if (sigmaSq <= 0)
        return -INFINITY;
    RealType x = point.first;
    RealType y = (alpha + 1.5) * std::log(sigmaSq);
    RealType degree = x - mu;
    degree *= degree;
    degree *= lambda;
    degree += 2 * beta;
    degree *= 0.5 / sigmaSq;
    y += degree;
    return pdfCoef - y;
}

template <typename RealType>
double NormalInverseGammaRand<RealType>::F(const Pair<RealType> &point) const
{
    RealType sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    RealType x = point.first;
    if (x == -INFINITY)
        return 0.0;
    double marginalY = this->Y.F(sigmaSq);
    if (x == INFINITY)
        return marginalY;
    if (sigmaSq == INFINITY)
        return this->X.F(x);
    if (x == mu)
        return 0.5 * marginalY;

    double lower = beta / sigmaSq;
    double normalCoef = (x - mu) * std::sqrt(lambda / beta);
    double logGammaAlpha = this->Y.GetLogGammaShape();
    double shape = alpha;
    double joint = math::integrate(
        [lower, normalCoef, logGammaAlpha, shape](double t) {
            double oneMinusT = 1.0 - t;
            double u = lower + t / oneMinusT;
            double logWeight = (shape - 1.0) * std::log(u) - u - logGammaAlpha;
            logWeight -= 2.0 * std::log(oneMinusT);
            double conditional = 0.5 * std::erfc(-normalCoef * std::sqrt(0.5 * u));
            return std::exp(logWeight) * conditional;
        },
        0.0, 1.0);
    return std::clamp(joint, 0.0, std::min(this->X.F(x), marginalY));
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
Pair<RealType> NormalInverseGammaRand<RealType>::operator()(G &g) const
{
    Pair<RealType> var;
    var.second = this->Y(g);
    double coef = std::sqrt(var.second / lambda);
    var.first = mu + coef * detail::standardNormal<RealType>(g);
    return var;
}

template <typename RealType>
long double NormalInverseGammaRand<RealType>::Correlation() const
{
    return (alpha > 2.0) ? 0.0 : NAN;
}

template <typename RealType>
Pair<RealType> NormalInverseGammaRand<RealType>::Mode() const
{
    return std::make_pair(mu, 2 * beta / (2 * alpha + 3));
}

} // namespace randlib

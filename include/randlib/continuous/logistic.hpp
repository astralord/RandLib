#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The LogisticRand class <BR>
 * Logistic distribution
 *
 * Notation: X ~ Logistic(μ, s)
 *
 * Related distributions: <BR>
 * 1 / (exp((X - μ) / s) + 1) ~ U(0, 1)
 */
template <typename RealType>
class LogisticRand final : public ContinuousDistribution<RealType>
{
    double mu = 0; ///< location μ
    double s = 1; ///< scale s
    double logS = 0; ///< log(s)

public:
    LogisticRand(double location = 0, double scale = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Infinite; }
    [[nodiscard]] RealType MinValue() const override { return -INFINITY; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    [[nodiscard]] inline double GetLocation() const noexcept { return mu; }
    [[nodiscard]] inline double GetScale() const noexcept { return s; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
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
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Entropy() const;

    /**
     * @fn FitLocation
     * fit location parameter via maximum-likelihood
     * @param sample
     */
    void FitLocation(std::span<const RealType> sample);
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template <typename RealType>
LogisticRand<RealType>::LogisticRand(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

template <typename RealType>
std::string LogisticRand<RealType>::Name() const
{
    return "Logistic(" + this->toStringWithPrecision(GetLocation()) + ", " + this->toStringWithPrecision(GetScale()) +
           ")";
}

template <typename RealType>
void LogisticRand<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Logistic distribution: location should be finite");
    mu = location;
}

template <typename RealType>
void LogisticRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Logistic distribution: scale of should be positive");
    s = scale;
    logS = std::log(s);
}

template <typename RealType>
double LogisticRand<RealType>::f(const RealType &x) const
{
    /// symmetric form exp(-|z|) / (1 + exp(-|z|))^2 never overflows
    double expMinusAbsZ = std::exp(-std::fabs(x - mu) / s);
    double denominator = 1 + expMinusAbsZ;
    return expMinusAbsZ / (denominator * denominator * s);
}

template <typename RealType>
double LogisticRand<RealType>::logf(const RealType &x) const
{
    double x0 = (mu - x) / s;
    double y = math::softplus(x0);
    y *= 2;
    y += logS;
    return x0 - y;
}

template <typename RealType>
double LogisticRand<RealType>::F(const RealType &x) const
{
    double expX = std::exp((mu - x) / s);
    return 1.0 / (1 + expX);
}

template <typename RealType>
double LogisticRand<RealType>::S(const RealType &x) const
{
    double expX = std::exp((mu - x) / s);
    return expX / (1 + expX);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType LogisticRand<RealType>::operator()(G &g) const
{
    /// there can be used rejection method from Laplace or Cauchy (Luc Devroye, p. 471) or ziggurat
    return mu + s * std::log(1.0 / detail::uniform01<RealType>(g) - 1);
}

template <typename RealType>
long double LogisticRand<RealType>::Mean() const
{
    return mu;
}

template <typename RealType>
long double LogisticRand<RealType>::Variance() const
{
    double sPi = s * randlib::consts::pi;
    return sPi * sPi / 3;
}

template <typename RealType>
std::complex<double> LogisticRand<RealType>::CFImpl(double t) const
{
    double pist = randlib::consts::pi * s * t;
    std::complex<double> y(0.0, t * mu);
    y = std::exp(y);
    y *= pist;
    y /= std::sinh(pist);
    return y;
}

template <typename RealType>
long double LogisticRand<RealType>::Entropy() const
{
    return 2 + logS;
}

template <typename RealType>
RealType LogisticRand<RealType>::quantileImpl(double p) const
{
    return mu - s * (std::log1pl(-p) - std::log(p));
}

template <typename RealType>
RealType LogisticRand<RealType>::quantileImpl1m(double p) const
{
    return mu - s * (std::log(p) - std::log1pl(-p));
}

template <typename RealType>
RealType LogisticRand<RealType>::Median() const
{
    return mu;
}

template <typename RealType>
RealType LogisticRand<RealType>::Mode() const
{
    return mu;
}

template <typename RealType>
long double LogisticRand<RealType>::Skewness() const
{
    return 0;
}

template <typename RealType>
long double LogisticRand<RealType>::ExcessKurtosis() const
{
    return 1.2;
}

template <typename RealType>
void LogisticRand<RealType>::FitLocation(std::span<const RealType> sample)
{
    double nHalf = 0.5 * sample.size();
    RealType root = 0;
    if (!math::findRootNewtonFirstOrder<RealType>(
            [this, sample, nHalf](RealType m) {
                double f1 = 0, f2 = 0;
                for (const RealType &x : sample) {
                    double aux = std::exp((m - x) / s);
                    double denom = 1.0 + aux;
                    f1 += 1.0 / denom;
                    denom *= denom;
                    f2 -= aux / denom;
                }
                f1 -= nHalf;
                return DoublePair(f1, f2);
            },
            root))
        throw this->fitError(FitErrorCode::Undefined, "Error in root-finding procedure");
    SetLocation(root);
}

} // namespace randlib

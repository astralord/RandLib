#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/gamma.hpp>

namespace randlib
{

/**
 * @brief The NakagamiDistribution class <BR>
 * Abstract class for Nakagami distribution
 *
 * Notation: X ~ Nakagami(μ, ω)
 *
 * Related distributions: <BR>
 * σX ~ Nakagami(μ, ωσ^2) <BR>
 * X^2 ~ Γ(μ, μ / ω)
 */
template <typename RealType>
class NakagamiDistribution : public ContinuousDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
    void sample(detail::EngineRef g, std::span<RealType> outputData) const override;

    double mu = 0.5; ///< shape μ
    double omega = 1; ///< spread ω
    GammaRand<RealType> Y{};
    double lgammaShapeRatio = 0; ///< log(Γ(μ + 0.5) / Γ(μ))

protected:
    NakagamiDistribution(double shape = 0.5, double spread = 1);

public:
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

protected:
    /**
     * @fn SetParameters
     * @param shape μ
     * @param spread ω
     */
    void SetParameters(double shape, double spread);

public:
    /**
     * @fn GetShape
     * @return shape μ
     */
    [[nodiscard]] inline double GetShape() const noexcept { return mu; }
    /**
     * @fn GetSpread
     * @return spread w
     */
    [[nodiscard]] inline double GetSpread() const noexcept { return omega; }
    /**
     * @fn GetLogGammaFunction
     * @return log(Γ(μ))
     */
    [[nodiscard]] inline double GetLogGammaFunction() const { return Y.GetLogGammaShape(); }
    /**
     * @fn GetLogGammaShapeRatio
     * @return log(Γ(μ + 0.5) / Γ(μ))
     */
    [[nodiscard]] inline double GetLogGammaShapeRatio() const noexcept { return lgammaShapeRatio; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double logf(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double FourthMoment() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

protected:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;
};

/**
 * @brief The NakagamiRand class <BR>
 * Nakagami distribution
 */
template <typename RealType>
class NakagamiRand final : public NakagamiDistribution<RealType>
{
public:
    NakagamiRand(double shape = 0.5, double spread = 1) : NakagamiDistribution<RealType>(shape, spread) {}
    [[nodiscard]] std::string Name() const override;
    using NakagamiDistribution<RealType>::SetParameters;
};

/**
 * @brief The ChiRand class <BR>
 * Chi distribution
 *
 * Notation: X ~ χ(k)
 *
 * Related distributions: <BR>
 * X ~ Nakagami(k/2, k) <BR>
 * X^2 ~ χ^2(k) <BR>
 * X^2 ~ Γ(k/2, 0.5)
 */
template <typename RealType>
class ChiRand final : public NakagamiDistribution<RealType>
{

public:
    explicit ChiRand(int degree);
    [[nodiscard]] std::string Name() const override;

    /**
     * @fn SetDegree
     * set degree k
     * @param degree
     */
    void SetDegree(int degree);
    /**
     * @fn GetDegree
     * @return degree k
     */
    [[nodiscard]] inline int GetDegree() const { return 2 * NakagamiDistribution<RealType>::GetShape(); }

    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
};

/**
 * @brief The MaxwellBoltzmannRand class <BR>
 * Maxwell-Boltzmann distribution
 *
 * Notation: X ~ MB(σ)
 *
 * Related distributions: <BR>
 * X / σ ~ χ(3) <BR>
 * X ~ Nakagami(1.5, 3σ^2)
 */
template <typename RealType>
class MaxwellBoltzmannRand final : public NakagamiDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

    double sigma = 1; ///< scale σ
public:
    explicit MaxwellBoltzmannRand(double scale);
    [[nodiscard]] std::string Name() const override;

    /**
     * @fn SetScale
     * set scale σ
     * @param scale
     */
    void SetScale(double scale);
    /**
     * @fn GetScale
     * @return scale σ
     */
    [[nodiscard]] double GetScale() const noexcept { return sigma; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;
};

/**
 * @brief The RayleighRand class <BR>
 * Rayleigh distribution
 *
 * Notation: X ~ Rayleigh(σ)
 *
 * Related distributions:
 * X / σ ~ χ(2)
 * X ~ Nakagami(1, 2σ^2)
 */
template <typename RealType>
class RayleighRand final : public NakagamiDistribution<RealType>
{
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }

    double sigma = 1; ///< scale σ
public:
    explicit RayleighRand(double scale = 1);
    [[nodiscard]] std::string Name() const override;

    /**
     * @fn SetScale
     * set scale σ
     * @param scale
     */
    void SetScale(double scale);
    /**
     * @fn GetScale
     * @return scale σ
     */
    [[nodiscard]] double GetScale() const noexcept { return sigma; }

    [[nodiscard]] double f(const RealType &x) const override;
    [[nodiscard]] double F(const RealType &x) const override;
    [[nodiscard]] double S(const RealType &x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

public:
    /**
     * @fn Fit
     * fit scale via maximum-likelihood method if unbiased = false,
     * otherwise set scale, returned by uniformly minimum variance unbiased estimator
     * @param sample
     */
    void Fit(std::span<const RealType> sample, bool unbiased = false);
};

} // namespace randlib

#include <randlib/continuous/exponential.hpp>
#include <randlib/continuous/normal.hpp>

namespace randlib
{

template <typename RealType>
NakagamiDistribution<RealType>::NakagamiDistribution(double shape, double spread)
{
    SetParameters(shape, spread);
}

template <typename RealType>
void NakagamiDistribution<RealType>::SetParameters(double shape, double spread)
{
    if (!(shape > 0.0))
        throw std::invalid_argument("Nakagami distribution: shape should be positive");
    if (!(spread > 0.0))
        throw std::invalid_argument("Nakagami distribution: spread should be positive");
    mu = shape;
    omega = spread;
    Y.SetParameters(mu, mu / omega);
    lgammaShapeRatio = std::lgammal(mu + 0.5) - Y.GetLogGammaShape();
}

template <typename RealType>
double NakagamiDistribution<RealType>::f(const RealType &x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0) {
        if (mu > 0.5)
            return 0.0;
        return (mu < 0.5) ? INFINITY : std::sqrt(randlib::consts::two_inv_pi / omega);
    }
    return 2 * x * Y.f(x * x);
}

template <typename RealType>
double NakagamiDistribution<RealType>::logf(const RealType &x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0) {
        if (mu > 0.5)
            return -INFINITY;
        return (mu < 0.5) ? INFINITY : 0.5 * (randlib::consts::ln2 - randlib::consts::lnpi - std::log(omega));
    }
    return randlib::consts::ln2 + std::log(x) + Y.logf(x * x);
}

template <typename RealType>
double NakagamiDistribution<RealType>::F(const RealType &x) const
{
    return (x > 0.0) ? Y.F(x * x) : 0.0;
}

template <typename RealType>
double NakagamiDistribution<RealType>::S(const RealType &x) const
{
    return (x > 0.0) ? Y.S(x * x) : 1.0;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType NakagamiDistribution<RealType>::operator()(G &g) const
{
    return std::sqrt(Y(g));
}

template <typename RealType>
void NakagamiDistribution<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    Y.Sample(g, outputData);
    for (RealType &var : outputData)
        var = std::sqrt(var);
}

template <typename RealType>
long double NakagamiDistribution<RealType>::Mean() const
{
    long double y = lgammaShapeRatio;
    y -= 0.5 * Y.GetLogRate();
    return std::exp(y);
}

template <typename RealType>
long double NakagamiDistribution<RealType>::Variance() const
{
    long double y = lgammaShapeRatio;
    y = std::exp(2 * y);
    return omega * (1 - y / mu);
}

template <typename RealType>
RealType NakagamiDistribution<RealType>::Median() const
{
    return std::sqrt(Y.Quantile(0.5));
}

template <typename RealType>
RealType NakagamiDistribution<RealType>::Mode() const
{
    long double mode = 1.0 - 0.5 / mu;
    if (mode <= 0.0)
        return 0.0;
    return std::sqrt(omega * mode);
}

template <typename RealType>
long double NakagamiDistribution<RealType>::Skewness() const
{
    long double thirdMoment = lgammaShapeRatio;
    thirdMoment -= 1.5 * Y.GetLogRate();
    thirdMoment = (mu + 0.5) * std::exp(thirdMoment);
    long double mean = Mean();
    long double variance = Variance();
    return (thirdMoment - mean * (3 * variance + mean * mean)) / std::pow(variance, 1.5);
}

template <typename RealType>
long double NakagamiDistribution<RealType>::FourthMoment() const
{
    long double fourthMoment = omega / mu;
    fourthMoment *= fourthMoment;
    fourthMoment *= mu * (mu + 1);
    return fourthMoment;
}

template <typename RealType>
long double NakagamiDistribution<RealType>::ExcessKurtosis() const
{
    long double mean = Mean();
    long double secondMoment = this->SecondMoment();
    long double thirdMoment = this->ThirdMoment();
    long double fourthMoment = FourthMoment();
    long double meanSq = mean * mean;
    long double variance = secondMoment - meanSq;
    long double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    long double denominator = variance * variance;
    return numerator / denominator - 3.0;
}

template <typename RealType>
RealType NakagamiDistribution<RealType>::quantileImpl(double p) const
{
    return std::sqrt(Y.Quantile(p));
}

template <typename RealType>
RealType NakagamiDistribution<RealType>::quantileImpl1m(double p) const
{
    return std::sqrt(Y.Quantile1m(p));
}

template <typename RealType>
std::complex<double> NakagamiDistribution<RealType>::CFImpl(double t) const
{
    if (mu >= 0.5)
        return ContinuousDistribution<RealType>::CFImpl(t);

    double re = this->ExpectedValue(
                    [t](double x) {
                        if (x == 0.0)
                            return 0.0;
                        return std::cos(t * x) - 1.0;
                    },
                    0, INFINITY) +
                1.0;

    double im = this->ExpectedValue([t](double x) { return std::sin(t * x); }, 0, INFINITY);

    return std::complex<double>(re, im);
}
/// NAKAGAMI
template <typename RealType>
std::string NakagamiRand<RealType>::Name() const
{
    return "Nakagami(" + this->toStringWithPrecision(this->GetShape()) + ", " +
           this->toStringWithPrecision(this->GetSpread()) + ")";
}
/// CHI
template <typename RealType>
ChiRand<RealType>::ChiRand(int degree)
{
    SetDegree(degree);
}

template <typename RealType>
std::string ChiRand<RealType>::Name() const
{
    return "Chi(" + this->toStringWithPrecision(GetDegree()) + ")";
}

template <typename RealType>
void ChiRand<RealType>::SetDegree(int degree)
{
    if (!(degree >= 1))
        throw std::invalid_argument("Chi distribution: degree parameter should be positive");
    NakagamiDistribution<RealType>::SetParameters(0.5 * degree, degree);
}

template <typename RealType>
long double ChiRand<RealType>::Skewness() const
{
    long double mean = this->Mean();
    long double sigmaSq = this->Variance();
    long double skew = mean * (1 - 2 * sigmaSq);
    skew /= std::pow(sigmaSq, 1.5);
    return skew;
}

template <typename RealType>
long double ChiRand<RealType>::ExcessKurtosis() const
{
    long double mean = this->Mean();
    long double sigmaSq = this->Variance();
    long double sigma = std::sqrt(sigmaSq);
    long double skew = Skewness();
    long double kurt = 1.0 - mean * sigma * skew;
    kurt /= sigmaSq;
    --kurt;
    return 2 * kurt;
}
/// MAXWELL-BOLTZMANN
template <typename RealType>
MaxwellBoltzmannRand<RealType>::MaxwellBoltzmannRand(double scale)
{
    SetScale(scale);
}

template <typename RealType>
std::string MaxwellBoltzmannRand<RealType>::Name() const
{
    return "Maxwell-Boltzmann(" + this->toStringWithPrecision(GetScale()) + ")";
}

template <typename RealType>
void MaxwellBoltzmannRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Maxwell-Boltzmann distribution: scale should be positive");
    sigma = scale;
    NakagamiDistribution<RealType>::SetParameters(1.5, 3 * sigma * sigma);
}

template <typename RealType>
double MaxwellBoltzmannRand<RealType>::f(const RealType &x) const
{
    if (x <= 0)
        return 0;
    double xAdj = x / sigma;
    double xAdjSq = xAdj * xAdj;
    double y = std::exp(-0.5 * xAdjSq);
    return randlib::consts::sqrt2 * randlib::consts::inv_sqrtpi * xAdjSq * y / sigma;
}

template <typename RealType>
double MaxwellBoltzmannRand<RealType>::F(const RealType &x) const
{
    if (x <= 0.0)
        return 0.0;
    double xAdj = randlib::consts::inv_sqrt2 * x / sigma;
    double y = std::exp(-xAdj * xAdj);
    y *= 2 * xAdj * randlib::consts::inv_sqrtpi;
    return std::erf(xAdj) - y;
}

template <typename RealType>
double MaxwellBoltzmannRand<RealType>::S(const RealType &x) const
{
    if (x <= 0.0)
        return 1.0;
    double xAdj = randlib::consts::inv_sqrt2 * x / sigma;
    double y = std::exp(-xAdj * xAdj);
    y *= 2 * xAdj * randlib::consts::inv_sqrtpi;
    return std::erfc(xAdj) + y;
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType MaxwellBoltzmannRand<RealType>::operator()(G &g) const
{
    RealType W = detail::standardExponential<RealType>(g);
    RealType N = detail::standardNormal<RealType>(g);
    return sigma * std::sqrt(2 * W + N * N);
}

template <typename RealType>
long double MaxwellBoltzmannRand<RealType>::Mean() const
{
    return 2 * randlib::consts::inv_sqrtpi * randlib::consts::sqrt2 * sigma;
}

template <typename RealType>
long double MaxwellBoltzmannRand<RealType>::Variance() const
{
    return (3 - 8.0 * randlib::consts::inv_pi) * sigma * sigma;
}

template <typename RealType>
RealType MaxwellBoltzmannRand<RealType>::Mode() const
{
    return randlib::consts::sqrt2 * sigma;
}

template <typename RealType>
long double MaxwellBoltzmannRand<RealType>::Skewness() const
{
    long double skewness = 3 * randlib::consts::pi - 8;
    skewness = 2.0 / skewness;
    skewness *= std::sqrt(skewness);
    return (16 - 5 * randlib::consts::pi) * skewness;
}

template <typename RealType>
long double MaxwellBoltzmannRand<RealType>::ExcessKurtosis() const
{
    long double numerator = 40 - 3 * randlib::consts::pi;
    numerator *= randlib::consts::pi;
    numerator -= 96;
    long double denominator = 3 * randlib::consts::pi - 8;
    denominator *= denominator;
    return 4 * numerator / denominator;
}
/// RAYLEIGH
template <typename RealType>
RayleighRand<RealType>::RayleighRand(double scale)
{
    SetScale(scale);
}

template <typename RealType>
std::string RayleighRand<RealType>::Name() const
{
    return "Rayleigh(" + this->toStringWithPrecision(GetScale()) + ")";
}

template <typename RealType>
void RayleighRand<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Rayleigh distribution: scale should be positive");
    sigma = scale;
    NakagamiDistribution<RealType>::SetParameters(1, 2 * sigma * sigma);
}

template <typename RealType>
double RayleighRand<RealType>::f(const RealType &x) const
{
    if (x <= 0)
        return 0.0;
    double y = x / (sigma * sigma);
    return y * std::exp(-0.5 * x * y);
}

template <typename RealType>
double RayleighRand<RealType>::F(const RealType &x) const
{
    if (x <= 0)
        return 0.0;
    double xAdj = x / sigma;
    return -std::expm1l(-0.5 * xAdj * xAdj);
}

template <typename RealType>
double RayleighRand<RealType>::S(const RealType &x) const
{
    if (x <= 0)
        return 1.0;
    double xAdj = x / sigma;
    return std::exp(-0.5 * xAdj * xAdj);
}

template <typename RealType>
template <std::uniform_random_bit_generator G>
RealType RayleighRand<RealType>::operator()(G &g) const
{
    double W = detail::standardExponential<RealType>(g);
    return sigma * std::sqrt(2 * W);
}

template <typename RealType>
long double RayleighRand<RealType>::Mean() const
{
    return sigma * randlib::consts::sqrtpi * randlib::consts::inv_sqrt2;
}

template <typename RealType>
long double RayleighRand<RealType>::Variance() const
{
    return (2.0 - randlib::consts::half_pi) * sigma * sigma;
}

template <typename RealType>
RealType RayleighRand<RealType>::quantileImpl(double p) const
{
    return sigma * std::sqrt(-2 * std::log1pl(-p));
}

template <typename RealType>
RealType RayleighRand<RealType>::quantileImpl1m(double p) const
{
    return sigma * std::sqrt(-2 * std::log(p));
}

template <typename RealType>
RealType RayleighRand<RealType>::Median() const
{
    static const double medianCoef = std::sqrt(randlib::consts::ln2 + randlib::consts::ln2);
    return sigma * medianCoef;
}

template <typename RealType>
RealType RayleighRand<RealType>::Mode() const
{
    return sigma;
}

template <typename RealType>
long double RayleighRand<RealType>::Skewness() const
{
    static const long double skewness =
        2 * randlib::consts::sqrtpi * (randlib::consts::pi - 3) / std::pow(4.0 - randlib::consts::pi, 1.5);
    return skewness;
}

template <typename RealType>
long double RayleighRand<RealType>::ExcessKurtosis() const
{
    static constexpr long double temp = 4 - randlib::consts::pi;
    static constexpr long double kurtosis =
        -(6 * randlib::consts::pi * randlib::consts::pi - 24 * randlib::consts::pi + 16) / (temp * temp);
    return kurtosis;
}

template <typename RealType>
void RayleighRand<RealType>::Fit(std::span<const RealType> sample, bool unbiased)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw this->fitError(FitErrorCode::WrongSample, this->POSITIVITY_VIOLATION);
    size_t n = sample.size();
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    double rawMoment = stats.second + stats.first * stats.first;
    double sigmaBiasedSq = 0.5 * rawMoment;
    if (unbiased == false) {
        SetScale(std::sqrt(sigmaBiasedSq));
    }
    /// Calculate unbiased sigma
    else if (n > 100) {
        double coef = 1.0 / (640 * std::pow(n, 5));
        coef -= 1.0 / (192 * std::pow(n, 3));
        coef += 0.125 / n;
        SetScale(std::exp(std::log1pl(coef) + 0.5 * sigmaBiasedSq)); /// err ~ o(n^{-6.5}) < 1e-13
    } else if (n > 15) {
        double coef = math::lfact(n);
        coef += math::lfact(n - 1);
        coef += 2 * n * randlib::consts::ln2;
        coef += 0.5 * std::log(n);
        coef -= math::lfact(2 * n);
        coef -= 0.5 * randlib::consts::lnpi;
        coef += 0.5 * std::log(sigmaBiasedSq);
        SetScale(std::exp(coef));
    } else {
        double scale = math::lfact(n - 1);
        scale += math::lfact(n);
        scale += 0.5 * std::log(n / randlib::consts::pi * sigmaBiasedSq);
        scale += 2 * n * randlib::consts::ln2;
        scale -= math::lfact(2 * n);
        SetScale(scale);
    }
}

} // namespace randlib

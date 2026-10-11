#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/exponential_family.hpp>
#include <randlib/continuous/geometric_stable.hpp>

namespace randlib
{

/**
 * @brief The AsymmetricLaplaceDistribution class <BR>
 * Abstract parent class for general Laplace and Asymmetric Laplace distributions
 */
template <typename RealType>
class AsymmetricLaplaceDistribution : public GeneralGeometricStableDistribution<RealType>
{
public:
    AsymmetricLaplaceDistribution(double shift = 0, double scale = 1, double asymmetry = 1);

protected:
    void ChangeLocation();

public:
    void SetScale(double scale);

    [[nodiscard]] inline double GetAsymmetry() const noexcept { return this->kappa; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
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

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Entropy() const;

    void FitScale(std::span<const RealType> sample);
};

/**
 * @brief The CenteredAsymmetricLaplaceRand class <BR>
 * Centered Laplace distribution
 *
 * Belongs to exponential family
 */
template <typename RealType>
class CenteredAsymmetricLaplaceRand final : public AsymmetricLaplaceDistribution<RealType>,
                                                           public ExponentialFamily<RealType, DoublePair>
{
public:
    CenteredAsymmetricLaplaceRand(double scale = 1, double asymmetry = 1) :
        AsymmetricLaplaceDistribution<RealType>(0.0, scale, asymmetry) {}
    [[nodiscard]] std::string Name() const override;

    using GeneralGeometricStableDistribution<RealType>::SetAsymmetry;

    /// Sufficient statistic (max(x, 0), max(-x, 0)) with natural parameters (-κ/γ, -1/(γκ))
    DoublePair SufficientStatistic(RealType x) const override;
    DoublePair SourceParameters() const override;
    DoublePair SourceToNatural(DoublePair sourceParameters) const override;
    double LogNormalizer(DoublePair theta) const override;
    DoublePair LogNormalizerGradient(DoublePair theta) const override;
    double CarrierMeasure(RealType) const override;
};

/**
 * @brief The ShiftedAsymmetricLaplaceDistribution class <BR>
 * Abstract parent class for shifted Laplace and Asymmetric Laplace distributions
 */
template <typename RealType>
class ShiftedAsymmetricLaplaceDistribution : public AsymmetricLaplaceDistribution<RealType>
{
public:
    ShiftedAsymmetricLaplaceDistribution(double shift = 0, double scale = 1, double asymmetry = 1) :
        AsymmetricLaplaceDistribution<RealType>(shift, scale, asymmetry) {}
    using GeneralGeometricStableDistribution<RealType>::SetShift;
    [[nodiscard]] inline double GetShift() const noexcept { return this->m; }

    void FitShift(std::span<const RealType> sample);
protected:
    void FitShiftAndScale(std::span<const RealType> sample);
};

/**
 * @brief The AsymmetricLaplaceRand class <BR>
 * Asymmetric Laplace distribution
 *
 * Notation: X ~ Asymmetric-Laplace(m, γ, κ)
 *
 * Related distributions: <BR>
 * X = m + γ * (Y / κ - W * κ), where Y, W ~ Exp(1) <BR>
 * X - m ~ GS(2, β, γ, γ(1 - κ^2) / κ) with arbitrary β
 */
template <typename RealType>
class AsymmetricLaplaceRand final : public ShiftedAsymmetricLaplaceDistribution<RealType>
{
public:
    AsymmetricLaplaceRand(double shift = 0, double scale = 1, double asymmetry = 1) :
        ShiftedAsymmetricLaplaceDistribution<RealType>(shift, scale, asymmetry) {}
    [[nodiscard]] std::string Name() const override;
    void SetAsymmetry(double asymmetry);
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(double asymmetry, G &g);

    using ShiftedAsymmetricLaplaceDistribution<RealType>::FitShiftAndScale;

private:
    DoublePair getOneSidedSums(std::span<const RealType> sample);
public:
    void FitAsymmetry(std::span<const RealType> sample);
};

/**
 * @brief The LaplaceRand class <BR>
 * Laplace distribution
 *
 * Notation: X ~ Laplace(m, γ)
 *
 * Related distributions: <BR>
 * X = m + γ * (Y - W), where Y, W ~ Exp(1) <BR>
 * X - m ~ GS(2, β, γ, 0) with arbitrary β
 */
template <typename RealType>
class LaplaceRand final : public ShiftedAsymmetricLaplaceDistribution<RealType>
{
public:
    LaplaceRand(double shift = 0, double scale = 1) :
        ShiftedAsymmetricLaplaceDistribution<RealType>(shift, scale, 1.0) {}
    [[nodiscard]] std::string Name() const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);
    void Fit(std::span<const RealType> sample) { this->FitShiftAndScale(sample); }
};

} // namespace randlib

#include <randlib/continuous/exponential.hpp>
#include <randlib/discrete/bernoulli.hpp>

namespace randlib
{

template < typename RealType >
AsymmetricLaplaceDistribution<RealType>::AsymmetricLaplaceDistribution(double shift, double scale, double asymmetry)
    : GeneralGeometricStableDistribution<RealType>(2.0, 0.0, scale, 0.0, shift)
{
    GeneralGeometricStableDistribution<RealType>::SetAsymmetry(asymmetry);
    ChangeLocation();
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::ChangeLocation()
{
    this->SetLocation((1.0 - this->kappaSq) * this->gamma * this->kappaInv);
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Laplace distribution: scale should be positive, but it's equal to "
                                    + std::to_string(scale));
    GeneralGeometricStableDistribution<RealType>::SetScale(scale);
    ChangeLocation();
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::f(const RealType & x) const
{
  return this->pdfLaplace(x - this->m);
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::logf(const RealType & x) const
{
  return this->logpdfLaplace(x - this->m);
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::F(const RealType & x) const
{
  return this->cdfLaplace(x - this->m);
}

template < typename RealType >
double AsymmetricLaplaceDistribution<RealType>::S(const RealType & x) const
{
  return this->cdfLaplaceCompl(x - this->m);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType AsymmetricLaplaceDistribution<RealType>::operator()(G &g) const
{
    RealType X = (this->kappa == 1) ? LaplaceRand<RealType>::StandardVariate(g)
                                    : AsymmetricLaplaceRand<RealType>::StandardVariate(this->kappa, g);
  return this->m + this->gamma * X;
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    if (this->kappa == 1) {
        for (RealType & var : outputData)
            var = this->m + this->gamma * LaplaceRand<RealType>::StandardVariate(g);
    }
    else {
        for (RealType & var : outputData)
            var = this->m + this->gamma * AsymmetricLaplaceRand<RealType>::StandardVariate(this->kappa, g);
    }
}

template < typename RealType >
std::complex<double> AsymmetricLaplaceDistribution<RealType>::CFImpl(double t) const
{
    double bt = this->gamma * t;
    double btSq = bt * bt;
    double denominator = (1 + this->kappaSq * btSq) * (1 + btSq / this->kappaSq);
    std::complex<double> y(std::cos(this->m * t), std::sin(this->m * t));
    std::complex<double> x(1, -this->kappa * bt), z(1, bt * this->kappaInv);
  return x * y * z / denominator;
}

template < typename RealType >
RealType AsymmetricLaplaceDistribution<RealType>::quantileImpl(double p) const
{
  return this->quantileLaplace(p);
}

template < typename RealType >
RealType AsymmetricLaplaceDistribution<RealType>::quantileImpl1m(double p) const
{
  return this->quantileLaplace1m(p);
}

template < typename RealType >
long double AsymmetricLaplaceDistribution<RealType>::Entropy() const
{
    double y = this->kappaInv + this->kappa;
  return std::log1pl(this->gamma * y);
}

template < typename RealType >
void AsymmetricLaplaceDistribution<RealType>::FitScale(std::span<const RealType> sample)
{
    double deviation = 0.0;
    for (const RealType & x : sample) {
        if (x > this->m)
            deviation += this->kappa * (x - this->m);
        else
            deviation -= (x - this->m) / this->kappa;
    }
    deviation /= sample.size();
    SetScale(deviation);
}

template < typename RealType >
void ShiftedAsymmetricLaplaceDistribution<RealType>::FitShift(std::span<const RealType> sample)
{
    /// Calculate median (considering asymmetry)
    /// we use root-finding algorithm for median search
    double minVar = *std::min_element(sample.begin(), sample.end());
    double maxVar = *std::max_element(sample.begin(), sample.end());
    double median = this->GetSampleMean(sample);

    if (!math::findRootNewtonFirstOrder<double>([this, sample] (double med)
    {
        double y = 0.0;
        for (const RealType & x : sample) {
            if (x > med)
                y -= this->kappaSq;
            else if (x < med)
                ++y;
        }
      return y;
    },
    minVar, maxVar, median
    ))
        throw this->fitError(FitErrorCode::Undefined, "Error in root-finding procedure");

    this->SetShift(median);
}

template < typename RealType >
void ShiftedAsymmetricLaplaceDistribution<RealType>::FitShiftAndScale(std::span<const RealType> sample)
{
    this->FitShift(sample);
    this->FitScale(sample);
}
template < typename RealType >
std::string CenteredAsymmetricLaplaceRand<RealType>::Name() const
{
  return "Centered Asymmetric-Laplace(" + this->toStringWithPrecision(this->GetScale()) + ", "
            + this->toStringWithPrecision(this->GetAsymmetry()) + ")";
}

template < typename RealType >
DoublePair CenteredAsymmetricLaplaceRand<RealType>::SufficientStatistic(RealType x) const
{
  return {std::max<double>(x, 0.0), std::max<double>(-x, 0.0)};
}

template < typename RealType >
DoublePair CenteredAsymmetricLaplaceRand<RealType>::SourceParameters() const
{
  return {this->GetScale(), this->GetAsymmetry()};
}

template < typename RealType >
DoublePair CenteredAsymmetricLaplaceRand<RealType>::SourceToNatural(DoublePair sourceParameters) const
{
    double scale = sourceParameters.first, asymmetry = sourceParameters.second;
  return {-asymmetry / scale, -1.0 / (scale * asymmetry)};
}

template < typename RealType >
double CenteredAsymmetricLaplaceRand<RealType>::LogNormalizer(DoublePair theta) const
{
  return std::log(-1.0 / theta.first - 1.0 / theta.second);
}

template < typename RealType >
DoublePair CenteredAsymmetricLaplaceRand<RealType>::LogNormalizerGradient(DoublePair theta) const
{
    double sum = -1.0 / theta.first - 1.0 / theta.second;
  return {1.0 / (theta.first * theta.first * sum), 1.0 / (theta.second * theta.second * sum)};
}

template < typename RealType >
double CenteredAsymmetricLaplaceRand<RealType>::CarrierMeasure(RealType) const
{
  return 0.0;
}
template < typename RealType >
std::string AsymmetricLaplaceRand<RealType>::Name() const
{
  return "Asymmetric-Laplace(" + this->toStringWithPrecision(this->GetShift()) + ", "
                                 + this->toStringWithPrecision(this->GetScale()) + ", "
                                 + this->toStringWithPrecision(this->GetAsymmetry()) + ")";
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::SetAsymmetry(double asymmetry)
{
    GeneralGeometricStableDistribution<RealType>::SetAsymmetry(asymmetry);
    this->ChangeLocation();
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType AsymmetricLaplaceRand<RealType>::StandardVariate(double asymmetry, G &g)
{
    RealType x = detail::standardExponential<RealType>(g) / asymmetry;
    RealType y = detail::standardExponential<RealType>(g) * asymmetry;
  return x - y;
}

template < typename RealType >
DoublePair AsymmetricLaplaceRand<RealType>::getOneSidedSums(std::span<const RealType> sample)
{
    double xPlus = 0.0, xMinus = 0.0;
    for (const RealType & x : sample) {
        if (x < this->m)
            xMinus += (x - this->m);
        else
            xPlus += (x - this->m);
    }
  return DoublePair(xPlus, xMinus);
}

template < typename RealType >
void AsymmetricLaplaceRand<RealType>::FitAsymmetry(std::span<const RealType> sample)
{
    auto [xPlus, xMinus] = getOneSidedSums(sample);
    double gammaN = this->gamma * sample.size();
    double asymmetry = 0;
    if (xPlus == xMinus)
        asymmetry = 1.0;
    else if (xPlus == 0)
        asymmetry = -xMinus / gammaN;
    else if (xMinus == 0) {
        asymmetry = gammaN / xPlus;
    }
    else {
        /// write down coefficients of quartic equation
        double a = xPlus / gammaN;
        double c = (xPlus + xMinus) / gammaN;
        double e = xMinus / gammaN;
        /// find coefficients for solution
        double ae = a * e;
        double delta1 = 2 * std::pow(c, 3) + 9 * c + 27 * (a + e) - 72 * ae * c;
        double delta0 = c * c + 3 + 12 * ae;
        double p2 = delta1 + std::sqrt(delta1 * delta1 - 4 * std::pow(delta0, 3));
        double Q = std::cbrt(0.5 * p2);
        double temp = 8 * a * a;
        double p = (8 * a * c - 3.0) / temp;
        double q = (-1.0 + 4 * a * c + temp) / (a * temp);
        double S = 0.5 * std::sqrt(-2.0 / 3 * p + (Q + delta0 / Q) / (3 * a));
        /// solve by general formula
        double c1 = -0.25 / a, b1 = -4 * S * S - 2 * p, b2 = q / S;
        if (b1 + b2 > 0)
            asymmetry = c1 + S + 0.5 * std::sqrt(b1 + b2);
        else
            asymmetry = c1 - S + 0.5 * std::sqrt(b1 - b2);
    }

    SetAsymmetry(asymmetry);
}
template < typename RealType >
std::string LaplaceRand<RealType>::Name() const
{
  return "Laplace(" + this->toStringWithPrecision(this->GetShift()) + ", "
                      + this->toStringWithPrecision(this->GetScale()) + ")";
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType LaplaceRand<RealType>::StandardVariate(G &g)
{
    RealType W = detail::standardExponential<RealType>(g);
  return BernoulliRand::StandardVariate(g) ? W : -W;
}

} // namespace randlib

#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>
#include <randlib/continuous/gamma.hpp>

namespace randlib
{

/**
 * @brief The InverseGammaRand class <BR>
 * Inverse-Gamma distribution
 *
 * X ~ Inv-Γ(α, β)
 *
 * Related distributions: <BR>
 * X = 1 / Y, where Y ~ Gamma(α, β)
 */
template <typename RealType>
class InverseGammaRand final : public ContinuousDistribution<RealType>
{
    double alpha = 1; ///< shape α
    double beta = 1; ///< rate β
    double pdfCoef = 0; ///< coefficient for faster pdf calculation

    GammaRand<RealType> X{};

public:
    InverseGammaRand(double shape = 1, double rate = 1);

    /// @return Inv-Γ(α, β) with shape α and rate (scale of the reciprocal's inverse) β
    [[nodiscard]] static InverseGammaRand FromShapeRate(double shape, double rate) { return InverseGammaRand(shape, rate); }

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::RightSemiFinite; }
    [[nodiscard]] RealType MinValue() const override { return 0; }
    [[nodiscard]] RealType MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double rate);
    [[nodiscard]] inline double GetShape() const noexcept { return alpha; }
    [[nodiscard]] inline double GetRate() const noexcept { return beta; }
    [[nodiscard]] inline double GetLogShape() const { return X.GetLogShape(); }
    [[nodiscard]] inline double GetLogRate() const { return X.GetLogRate(); }

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
    [[nodiscard]] double GetLogGammaShape() const { return X.GetLogGammaShape(); }
};

} // namespace randlib

namespace randlib
{

template < typename RealType >
InverseGammaRand<RealType>::InverseGammaRand(double shape, double rate)
{
    SetParameters(shape, rate);
}

template < typename RealType >
std::string InverseGammaRand<RealType>::Name() const
{
  return "Inverse-Gamma(" + this->toStringWithPrecision(GetShape()) + ", " + this->toStringWithPrecision(GetRate()) + ")";
}

template < typename RealType >
void InverseGammaRand<RealType>::SetParameters(double shape, double rate)
{
    if (!(shape > 0.0))
        throw std::invalid_argument("Inverse-Gamma distribution: shape should be positive");
    if (!(rate > 0.0))
        throw std::invalid_argument("Inverse-Gamma distribution: rate should be positive");
    X.SetParameters(shape, rate);
    alpha = X.GetShape();
    beta = X.GetRate();
    pdfCoef = alpha * X.GetLogRate() - X.GetLogGammaShape();
}

template < typename RealType >
double InverseGammaRand<RealType>::logf(const RealType & x) const
{
    if (x <= 0.0)
      return -INFINITY;
    double logX = std::log(x);
    double y = -(alpha - 1.0) * logX;
    y -= beta / x;
    y += pdfCoef;
  return y - 2 * logX;
}

template < typename RealType >
double InverseGammaRand<RealType>::F(const RealType & x) const
{
  return (x > 0.0) ? X.S(1.0 / x) : 0.0;
}

template < typename RealType >
double InverseGammaRand<RealType>::S(const RealType & x) const
{
  return (x > 0.0) ? X.F(1.0 / x) : 1.0;
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType InverseGammaRand<RealType>::operator()(G &g) const
{
  return 1.0 / X(g);
}

template < typename RealType >
void InverseGammaRand<RealType>::sample(detail::EngineRef g, std::span<RealType> outputData) const
{
    X.Sample(g, outputData);
    for (RealType &var : outputData)
        var = 1.0 / var;
}



template < typename RealType >
long double InverseGammaRand<RealType>::Mean() const
{
  return (alpha > 1) ? beta / (alpha - 1) : INFINITY;
}

template < typename RealType >
long double InverseGammaRand<RealType>::Variance() const
{
    if (alpha <= 2)
      return INFINITY;
    double var = beta / (alpha - 1);
    var *= var;
  return var / (alpha - 2);
}

template < typename RealType >
RealType InverseGammaRand<RealType>::Median() const
{
  return 1.0 / X.Median();
}

template < typename RealType >
RealType InverseGammaRand<RealType>::Mode() const
{
  return beta / (alpha + 1);
}

template < typename RealType >
RealType InverseGammaRand<RealType>::quantileImpl(double p) const
{
  return 1.0 / X.Quantile1m(p);
}

template < typename RealType >
RealType InverseGammaRand<RealType>::quantileImpl1m(double p) const
{
  return 1.0 / X.Quantile(p);
}

template < typename RealType >
long double InverseGammaRand<RealType>::Skewness() const
{
  return (alpha > 3) ? 4 * std::sqrt(alpha - 2) / (alpha - 3) : INFINITY;
}

template < typename RealType >
long double InverseGammaRand<RealType>::ExcessKurtosis() const
{
    if (alpha <= 4)
      return INFINITY;
    long double numerator = 30 * alpha - 66.0;
    long double denominator = (alpha - 3) * (alpha - 4);
  return numerator / denominator;
}

} // namespace randlib

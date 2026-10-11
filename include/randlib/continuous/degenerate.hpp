#pragma once

#include <span>

#include <algorithm>
#include <functional>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The DegenerateRand class <BR>
 * Degenerate distribution
 *
 * f(x|a) = δ(a)
 *
 * Notation: X ~ δ(a)
 */
class DegenerateRand final : public ContinuousDistribution<double>
{
    double a = 0; ///< value

public:
    explicit DegenerateRand(double value = 0);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] double MinValue() const override { return a; }
    [[nodiscard]] double MaxValue() const override { return a; }

    void SetValue(double value);
    [[nodiscard]] inline double GetValue() const noexcept { return a; }

    [[nodiscard]] double f(const double & x) const override;
    [[nodiscard]] double logf(const double & x) const override;
    [[nodiscard]] double F(const double & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] double operator()(G &g) const;
protected:
    double variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] double Median() const override;
    [[nodiscard]] double Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] long double Entropy() const;

    /**
     * @fn Fit
     * @param sample
     */
    void Fit(std::span<const double> sample);
};

} // namespace randlib

namespace randlib
{

inline DegenerateRand::DegenerateRand(double value)
{
    SetValue(value);
}

inline std::string DegenerateRand::Name() const
{
  return "Degenerate";
}

inline void DegenerateRand::SetValue(double value)
{
    if (!std::isfinite(value))
        throw std::invalid_argument("Degenerate distribution: value should be finite");
    a = value;
}

inline double DegenerateRand::f(const double & x) const
{
  return (x == a) ? INFINITY : 0.0;
}

inline double DegenerateRand::logf(const double & x) const
{
  return (x == a) ? INFINITY : -INFINITY;
}

inline double DegenerateRand::F(const double & x) const
{
  return (x < a) ? 0.0 : 1.0;
}

template <std::uniform_random_bit_generator G>
double DegenerateRand::operator()(G &) const
{
  return a;
}

inline long double DegenerateRand::Mean() const
{
  return a;
}

inline long double DegenerateRand::Variance() const
{
  return 0.0;
}

inline std::complex<double> DegenerateRand::CFImpl(double t) const
{
    double re = std::cos(a * t);
    double im = std::sin(a * t);
  return std::complex<double>(re, im);
}

inline double DegenerateRand::quantileImpl(double) const
{
  return a;
}

inline double DegenerateRand::quantileImpl1m(double) const
{
  return a;
}

inline double DegenerateRand::Median() const
{
  return a;
}

inline double DegenerateRand::Mode() const
{
  return a;
}

inline long double DegenerateRand::Skewness() const
{
  return NAN;
}

inline long double DegenerateRand::ExcessKurtosis() const
{
  return NAN;
}

inline long double DegenerateRand::Entropy() const
{
  return 0.0;
}

inline void DegenerateRand::Fit(std::span<const double> sample)
{
    if (sample.empty())
        throw this->fitError(FitErrorCode::TooFewElements, "The sample is empty");
    if (std::adjacent_find(sample.begin(), sample.end(), std::not_equal_to<>()) != sample.end())
        throw this->fitError(FitErrorCode::WrongSample, "All elements should be equal to each other");
    SetValue(sample.front());
}

} // namespace randlib

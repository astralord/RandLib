#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The RademacherRand class <BR>
 * Rademacher distribution
 *
 * P(X = k) = 0.5 * 1_{|k| = 1}
 *
 * Notation: X ~ Rademacher
 *
 * Related distributions: <BR>
 * If Y ~ Bernoulli(0.5), then 2Y - 1 ~ Rademacher
 */
class RademacherRand final : public DiscreteDistribution<int>
{
public:
    RademacherRand();
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] int MinValue() const override { return -1; }
    [[nodiscard]] int MaxValue() const override { return 1; }

    [[nodiscard]] double P(const int & k) const override;
    [[nodiscard]] double logP(const int & k) const override;
    [[nodiscard]] double F(const int & k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] int operator()(G &g) const;
protected:
    int variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] int Median() const override;
    [[nodiscard]] int Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    int quantileImpl(double p) const override;
    int quantileImpl1m(double p) const override;
    std::complex<double> CFImpl(double t) const override;

public:
    [[nodiscard]] double Entropy();
};

} // namespace randlib

#include <randlib/discrete/bernoulli.hpp>

namespace randlib
{

inline RademacherRand::RademacherRand()
{
}

inline std::string RademacherRand::Name() const
{
  return "Rademacher";
}

inline double RademacherRand::P(const int & k) const
{
  return (k == 1 || k == -1) ? 0.5 : 0.0;
}

inline double RademacherRand::logP(const int & k) const
{
  return (k == 1 || k == -1) ? -randlib::consts::ln2 : -INFINITY;
}

inline double RademacherRand::F(const int & k) const
{
    if (k < -1)
      return 0;
  return (k < 1) ? 0.5 : 1.0;
}

template <std::uniform_random_bit_generator G>
int RademacherRand::operator()(G &g) const
{
  return BernoulliRand::StandardVariate(g) ? 1 : -1;
}

inline long double RademacherRand::Mean() const
{
  return 0;
}

inline long double RademacherRand::Variance() const
{
  return 1;
}

inline std::complex<double> RademacherRand::CFImpl(double t) const
{
  return std::cos(t);
}

inline int RademacherRand::Median() const
{
  return -1;
}

inline int RademacherRand::Mode() const
{
    /// both -1 and 1 are modes; the smaller one, like the median
  return -1;
}

inline long double RademacherRand::Skewness() const
{
  return 0.0;
}

inline long double RademacherRand::ExcessKurtosis() const
{
  return -2.0;
}

inline int RademacherRand::quantileImpl(double p) const
{
  return (p <= 0.5) ? -1 : 1;
}

inline int RademacherRand::quantileImpl1m(double p) const
{
  return (p >= 0.5) ? -1 : 1;
}

inline double RademacherRand::Entropy()
{
  return randlib::consts::ln2;
}

} // namespace randlib

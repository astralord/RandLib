#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/univariate.hpp>

namespace randlib
{

/**
 * @brief The SingularDistribution class <BR>
 * Abstract class for all singular distributions
 */
class SingularDistribution : public UnivariateDistribution<double>
{
protected:
    SingularDistribution();
public:
    virtual ~SingularDistribution() {}
protected:

private:
    [[nodiscard]] double Hazard(const double&) const override;
    [[nodiscard]] double Mode() const override;
    long double ExpectedValue(const std::function<double (double)> &funPtr, double minPoint, double maxPoint) const override;
    [[nodiscard]] double LikelihoodFunction(std::span<const double> sample) const override;
    [[nodiscard]] double LogLikelihoodFunction(std::span<const double> sample) const override;
};

} // namespace randlib

namespace randlib
{

inline SingularDistribution::SingularDistribution()
{

}

inline double SingularDistribution::Hazard(const double &) const
{
  return NAN;
}

inline long double SingularDistribution::ExpectedValue(const std::function<double (double)> &, double, double ) const
{
  return NAN;
}

inline double SingularDistribution::Mode() const
{
  return NAN;
}

inline double SingularDistribution::LikelihoodFunction(std::span<const double>) const
{
  return NAN;
}

inline double SingularDistribution::LogLikelihoodFunction(std::span<const double>) const
{
  return NAN;
}

} // namespace randlib

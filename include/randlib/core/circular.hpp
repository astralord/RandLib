#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/continuous.hpp>

namespace randlib
{

/**
 * @brief The CircularDistribution class <BR>
 * Abstract class for all continuous circular distributions
 *
 * Circular summaries are represented by circular moments rather than linear moments.
 */
template <typename RealType>
class CircularDistribution : public ContinuousDistribution<RealType>
{
protected:
    double loc{};

    CircularDistribution(double location = 0);

public:
    virtual ~CircularDistribution() {}

protected:
public:
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return loc - randlib::consts::pi; }
    [[nodiscard]] RealType MaxValue() const override { return loc + randlib::consts::pi; }

    void SetLocation(double location);
    [[nodiscard]] inline double GetLocation() const noexcept { return loc; }

    [[nodiscard]] long double Mean() const override { return NAN; }
    [[nodiscard]] long double Variance() const override { return NAN; }
    [[nodiscard]] long double Skewness() const override { return NAN; }
    [[nodiscard]] long double ExcessKurtosis() const override { return NAN; }

    /**
     * @return Circular mean of random variable
     */
    [[nodiscard]] virtual long double CircularMean() const = 0;

    /**
     * @return Circular variance of random variable
     */
    [[nodiscard]] virtual long double CircularVariance() const = 0;
};

} // namespace randlib

namespace randlib
{

template <typename RealType>
CircularDistribution<RealType>::CircularDistribution(double location)
{
    SetLocation(location);
}

template <typename RealType>
void CircularDistribution<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Circular distribution: location should be finite");
    loc = location;
}

} // namespace randlib

#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/bivariate.hpp>

namespace randlib
{

/**
 * @brief The ContinuousBivariateDistribution class <BR>
 * Abstract class for all bivariate probability distributions
 * with marginal continuous distributions
 */
template <class T1, class T2, typename RealType>
class ContinuousBivariateDistribution : public BivariateDistribution< T1, T2, RealType>
{
    static_assert(std::is_base_of_v<ContinuousDistribution<RealType>, T1>, "T1 must be a descendant of ContinuousDistribution");
    static_assert(std::is_base_of_v<ContinuousDistribution<RealType>, T2>, "T2 must be a descendant of ContinuousDistribution");

protected:
    ContinuousBivariateDistribution() {}
public:
    virtual ~ContinuousBivariateDistribution() {}
protected:

public:
    [[nodiscard]] virtual double f(const Pair<RealType> &point) const = 0;
    [[nodiscard]] virtual double logf(const Pair<RealType> &point) const = 0;
};

} // namespace randlib

#include <randlib/core/continuous.hpp>
#include <randlib/continuous/student_t.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/inverse_gamma.hpp>

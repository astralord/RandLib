#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/bivariate.hpp>

namespace randlib
{

/**
 * @brief The DiscreteBivariateDistribution class <BR>
 * Abstract class for all bivariate probability distributions
 * with marginal discrete distributions
 */
template <class T1, class T2, typename IntType>
class DiscreteBivariateDistribution : public BivariateDistribution<T1, T2, IntType>
{
    static_assert(std::is_base_of_v<DiscreteDistribution<IntType>, T1>,
                  "T1 must be a descendant of DiscreteDistribution");
    static_assert(std::is_base_of_v<DiscreteDistribution<IntType>, T2>,
                  "T2 must be a descendant of DiscreteDistribution");

protected:
    DiscreteBivariateDistribution() {}

public:
    virtual ~DiscreteBivariateDistribution() {}

protected:
public:
    [[nodiscard]] virtual double P(const Pair<IntType> &point) const = 0;
    [[nodiscard]] virtual double logP(const Pair<IntType> &point) const = 0;
};

} // namespace randlib

#include <randlib/core/discrete.hpp>
#include <randlib/discrete/binomial.hpp>

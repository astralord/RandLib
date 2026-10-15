#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/distribution.hpp>

namespace randlib
{

/**
 * @brief The BivariateDistribution class <BR>
 * Abstract class for all bivariate probability distributions
 */
template <class T1, class T2, typename T>
class BivariateDistribution : public ProbabilityDistribution<Pair<T>>
{
protected:
    T1 X{}; ///< 1st marginal distributions
    T2 Y{}; ///< 2nd marginal distributions
    BivariateDistribution() {}

public:
    virtual ~BivariateDistribution() {}

protected:
public:
    [[nodiscard]] Pair<T> MinValue() const override { return Pair<T>(X.MinValue(), Y.MinValue()); }
    [[nodiscard]] Pair<T> MaxValue() const override { return Pair<T>(X.MaxValue(), Y.MaxValue()); }

    [[nodiscard]] virtual LongDoublePair Mean() const final;
    virtual LongDoubleTriplet Covariance() const final;
    virtual long double Correlation() const = 0;
    [[nodiscard]] virtual std::pair<T1, T2> GetMarginalDistributions() const final;
    [[nodiscard]] virtual Pair<T> Mode() const = 0;
};

} // namespace randlib

#include <randlib/core/univariate.hpp>
#include <randlib/discrete/binomial.hpp>
#include <randlib/continuous/student_t.hpp>
#include <randlib/continuous/normal.hpp>
#include <randlib/continuous/inverse_gamma.hpp>

namespace randlib
{

template <class T1, class T2, typename T>
LongDoublePair BivariateDistribution<T1, T2, T>::Mean() const
{
    return std::make_pair(X.Mean(), Y.Mean());
}

template <class T1, class T2, typename T>
LongDoubleTriplet BivariateDistribution<T1, T2, T>::Covariance() const
{
    double var1 = X.Variance();
    double var2 = Y.Variance();
    double corr = Correlation() * var1 * var2;
    return std::make_tuple(var1, corr, var2);
}

template <class T1, class T2, typename T>
std::pair<T1, T2> BivariateDistribution<T1, T2, T>::GetMarginalDistributions() const
{
    return std::make_pair(X, Y);
}

} // namespace randlib

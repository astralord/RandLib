#ifndef DISCRETEBIVARIATEDISTRIBUTION_H
#define DISCRETEBIVARIATEDISTRIBUTION_H

#include "BivariateDistribution.h"
#include "../univariate/discrete/DiscreteDistribution.h"

/**
 * @brief The DiscreteBivariateDistribution class <BR>
 * Abstract class for all bivariate probability distributions
 * with marginal discrete distributions
 */
template < class T1, class T2, typename IntType >
class RANDLIBSHARED_EXPORT DiscreteBivariateDistribution : public BivariateDistribution< T1, T2, IntType >
{
    static_assert(std::is_base_of_v<DiscreteDistribution<IntType>, T1>, "T1 must be a descendant of DiscreteDistribution");
    static_assert(std::is_base_of_v<DiscreteDistribution<IntType>, T2>, "T2 must be a descendant of DiscreteDistribution");

protected:
    DiscreteBivariateDistribution() {}
    virtual ~DiscreteBivariateDistribution() {}

public:
    virtual double P(const Pair<IntType> &point) const = 0;
    virtual double logP(const Pair<IntType> &point) const = 0;
};


#endif // DISCRETEBIVARIATEDISTRIBUTION_H

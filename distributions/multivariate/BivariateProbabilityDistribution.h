#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"
#include "../univariate/UnivariateProbabilityDistribution.h"
#include "../../math/SquareMatrix.h"

/**
 * @brief The BivariateProbabilityDistribution class
 */
template < typename T1, typename T2 >
class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution< std::pair<T1, T2> >
{
    typedef std::pair<T1, T2> Pair;

public:
    BivariateProbabilityDistribution();
    virtual ~BivariateProbabilityDistribution() {}

    virtual double f(Pair point) const = 0;

    virtual DoublePair Mean() const = 0;
    virtual DoubleTriplet Covariance() const = 0;
    virtual double Correlation() const = 0;

    virtual void GetMarginalDistributions(UnivariateProbabilityDistribution<T1> &distribution1, UnivariateProbabilityDistribution<T2> &distribution2) const = 0;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

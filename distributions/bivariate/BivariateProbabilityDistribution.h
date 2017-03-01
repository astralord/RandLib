#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"
#include "../univariate/continuous/ContinuousDistribution.h"

/**
 * @brief The BivariateProbabilityDistribution class
 */
class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution< DoublePair >
{

public:
    BivariateProbabilityDistribution() {}
    virtual ~BivariateProbabilityDistribution() {}
    virtual double f(const DoublePair &point) const = 0;

    virtual DoublePair Mean() const = 0;
    virtual DoubleTriplet Covariance() const = 0;
    virtual double Correlation() const = 0;

    virtual void GetMarginalDistributions(ContinuousDistribution &distribution1, ContinuousDistribution &distribution2) const = 0;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

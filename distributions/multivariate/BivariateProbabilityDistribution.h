#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"
#include "../../math/Matrix.h"

/**
 * @brief The BivariateProbabilityDistribution class
 */
class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution<DoublePair>
{
public:
    BivariateProbabilityDistribution();
    virtual ~BivariateProbabilityDistribution() {}

    virtual double f(DoublePair point) const = 0;

    virtual bool Covariance(Matrix &matrix) const = 0;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

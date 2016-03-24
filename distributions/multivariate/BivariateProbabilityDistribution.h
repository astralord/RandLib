#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"
#include "../../math/Matrix.h"

/**
 * @brief The BivariateProbabilityDistribution class
 */
class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution<double2d>
{
public:
    BivariateProbabilityDistribution();
    virtual ~BivariateProbabilityDistribution() {}

    virtual double f(double2d point) const = 0;

    virtual bool Covariance(Matrix &matrix) const = 0;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

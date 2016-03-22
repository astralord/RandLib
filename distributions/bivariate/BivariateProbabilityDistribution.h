#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"

/**
 * @brief The BivariateProbabilityDistribution class
 */
class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution<double2d>
{
public:
    BivariateProbabilityDistribution();
    virtual ~BivariateProbabilityDistribution() {}

    virtual double f(double2d grid) const = 0;

    // TODO: add virtual Matrix Covariance() const = 0;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

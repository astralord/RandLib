#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"

class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution<double2d>
{
public:
    BivariateProbabilityDistribution();

    /**
     * @brief Covariance
     * @return Covariance of random variable
     */
    virtual double2d Covariance() const = 0;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

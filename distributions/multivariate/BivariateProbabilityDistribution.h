#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"
#include "../univariate/UnivariateProbabilityDistribution.h"
#include "../../math/SquareMatrix.h"

/**
 * @brief The BivariateProbabilityDistribution class
 */
class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution<DoublePair>
{
public:
    BivariateProbabilityDistribution();
    virtual ~BivariateProbabilityDistribution() {}

    virtual double f(DoublePair point) const = 0;

    virtual void Covariance(SquareMatrix<2> &matrix) const = 0;

    virtual double Correlation() const = 0;

    virtual void getFirstMarginalDistribution(UnivariateProbabilityDistribution &distribution) const = 0;
    virtual void getSecondMarginalDistribution(UnivariateProbabilityDistribution &distribution) const = 0;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

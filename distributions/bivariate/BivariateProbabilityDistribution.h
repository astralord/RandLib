#ifndef BIVARIATEPROBABILITYDISTRIBUTION_H
#define BIVARIATEPROBABILITYDISTRIBUTION_H

#include "../ProbabilityDistribution.h"
#include "../univariate/continuous/ContinuousDistribution.h"

/**
 * @brief The BivariateProbabilityDistribution class <BR>
 * Abstract class for all bivariate probability distributions
 */
template < class T1, class T2 >
class RANDLIBSHARED_EXPORT BivariateProbabilityDistribution : public ProbabilityDistribution< DoublePair >
{
protected:
    T1 X{}; ///< 1st marginal distributions
    T2 Y{}; ///< 2nd marginal distributions
    BivariateProbabilityDistribution() {}
    virtual ~BivariateProbabilityDistribution() {}

public:
    virtual double f(const DoublePair &point) const = 0;
    virtual double logf(const DoublePair &point) const = 0;

    virtual DoublePair Mean() const final;
    virtual DoubleTriplet Covariance() const final;
    virtual double Correlation() const = 0;

    virtual std::pair<T1, T2> GetMarginalDistributions() const final;
};

#endif // BIVARIATEPROBABILITYDISTRIBUTION_H

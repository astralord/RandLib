#ifndef BIVARIATEDISTRIBUTION_H
#define BIVARIATEDISTRIBUTION_H

#include "../ProbabilityDistribution.h"
#include "../univariate/UnivariateDistribution.h"

/**
 * @brief The BivariateDistribution class <BR>
 * Abstract class for all bivariate probability distributions
 */
template < class T1, class T2, typename T >
class RANDLIBSHARED_EXPORT BivariateDistribution : public ProbabilityDistribution< Pair<T> >
{
protected:
    T1 X{}; ///< 1st marginal distributions
    T2 Y{}; ///< 2nd marginal distributions
    BivariateDistribution() {}
    virtual ~BivariateDistribution() {}

public:
    Pair<T> MinValue() const { return Pair<T>(X.MinValue(), Y.MinValue()); }
    Pair<T> MaxValue() const { return Pair<T>(X.MaxValue(), Y.MaxValue()); }

    void Reseed(unsigned long seed) const override;

    virtual LongDoublePair Mean() const final;
    virtual LongDoubleTriplet Covariance() const final;
    virtual long double Correlation() const = 0;
    virtual std::pair<T1, T2> GetMarginalDistributions() const final;
    virtual Pair<T> Mode() const = 0;
};

#endif // BIVARIATEDISTRIBUTION_H

#ifndef COMPOUNDCONTINUOUSRAND_H
#define COMPOUNDCONTINUOUSRAND_H

#include "ContinuousDistribution.h"
#include "../CompoundDistribution.h"

/**
 * @brief The CompoundContinuousRand class
 */
template <typename T>
class RANDLIBSHARED_EXPORT CompoundContinuousRand : public CompoundDistribution<double, T>, ContinuousDistribution
{
public:
    CompoundContinuousRand(const ContinuousDistribution &X, const ContinuousDistribution &Y) :
        CompoundContinuousRand(X, Y)
    {}
    virtual ~CompoundContinuousRand() {}

    double f(double x) const override;
};

#endif // COMPOUNDCONTINUOUSRAND_H

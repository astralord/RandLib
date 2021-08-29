#ifndef EXPONENTIALFAMILY_H
#define EXPONENTIALFAMILY_H

#include "UnivariateDistribution.h"

template < typename T, typename P >
class RANDLIBSHARED_EXPORT ExponentialFamily
{
public:
    ExponentialFamily() {}
    virtual ~ExponentialFamily() {}

    virtual P SufficientStatistic(T x) const = 0;
    virtual P SourceParameters() const = 0;
    virtual P SourceToNatural(P sourceParameters) const = 0;
    virtual P NaturalParameters() const;

    virtual double LogNormalizer(P theta) const = 0;
    virtual P LogNormalizerGradient(P theta) const = 0;
    virtual double CarrierMeasure(T x) const = 0;

    double ProbabilityMeasure(T x) const;
    double LogProbabilityMeasure(T x) const;

    double KullbackLeiblerDivergence(P parameters) const;
    virtual double CrossEntropyAdjusted(P parameters) const;
    virtual double EntropyAdjusted() const;
};

#endif // EXPONENTIALFAMILY_H

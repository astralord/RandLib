#ifndef COMPOUNDRAND_H
#define COMPOUNDRAND_H

#include "UnivariateProbabilityDistribution.h"

/**
 * @brief The CompoundRand class
 */
template < typename T1, typename T2 >
class RANDLIBSHARED_EXPORT CompoundRand : public UnivariateProbabilityDistribution<T1>
{
    const UnivariateProbabilityDistribution<T1> &X;
    const UnivariateProbabilityDistribution<T2> &Y;
public:
    CompoundRand(const UnivariateProbabilityDistribution<T1> & leftRV, const UnivariateProbabilityDistribution<T2> & rightRV);
    virtual ~CompoundRand() {}

    T1 variate() const override;

    double Mean() const override;
    double Variance() const override;
    std::complex<double> CF(double t) const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // COMPOUNDRAND_H

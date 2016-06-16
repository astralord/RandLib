#ifndef COMPOUNDDISTRIBUTION_H
#define COMPOUNDDISTRIBUTION_H

#include "UnivariateProbabilityDistribution.h"

/**
 * @brief The CompoundDistribution class
 */
template < typename T1, typename T2 >
class RANDLIBSHARED_EXPORT CompoundDistribution : public UnivariateProbabilityDistribution<T1>
{
protected:
    const UnivariateProbabilityDistribution<T1> &X;
    const UnivariateProbabilityDistribution<T2> &Y;
public:
    CompoundDistribution(const UnivariateProbabilityDistribution<T1> & leftRV, const UnivariateProbabilityDistribution<T2> & rightRV);
    virtual ~CompoundDistribution() {}
    std::string name() const override;

    double F(T1 x) const override;
    T1 variate() const override;

    double Mean() const override;
    double Variance() const override;
    std::complex<double> CF(double t) const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // COMPOUNDDISTRIBUTION_H

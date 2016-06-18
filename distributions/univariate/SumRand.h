#ifndef SUMRAND_H
#define SUMRAND_H

#include "UnivariateProbabilityDistribution.h"
#include "continuous/ContinuousDistribution.h"
#include "discrete/DiscreteDistribution.h"

/**
 * @brief The SumContinuousRand class
 */
template < typename T>
class RANDLIBSHARED_EXPORT SumContinuousRand : public ContinuousDistribution
{
    const ContinuousDistribution &X;
    const UnivariateProbabilityDistribution<T> &Y;
public:
    SumContinuousRand(const ContinuousDistribution & leftRV, const UnivariateProbabilityDistribution<T> & rightRV);
    virtual ~SumContinuousRand() {}
    std::string name() const override;

    SUPPORT_TYPE supportType() const override;
    double MinValue() const override;
    double MaxValue() const override;

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;
    std::complex<double> CF(double t) const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

const SumContinuousRand<double> operator+(const ContinuousDistribution& left, const ContinuousDistribution& right) {
    return SumContinuousRand<double>(left, right);
}


#endif // SUMRAND_H

#ifndef CIRCULARDISTRIBUTION_H
#define CIRCULARDISTRIBUTION_H

#include "../ContinuousDistribution.h"

/**
 * @brief The CircularDistribution class <BR>
 * Abstract class for all continuous circular distributions
 *
 * Note that all the moments are now useless, we implement circular moments instead.
 */
template < typename RealType = double >
class RANDLIBSHARED_EXPORT CircularDistribution : public ContinuousDistribution<RealType>
{
protected:
    double loc{};

    CircularDistribution(double location = 0);
    virtual ~CircularDistribution() {}

public:
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    RealType MinValue() const override { return loc - M_PI; }
    RealType MaxValue() const override { return loc + M_PI; }

    void SetLocation(double location);
    inline double GetLocation() const { return loc; }

    long double Mean() const override { return NAN; }
    long double Variance() const override { return NAN; }
    long double Skewness() const override { return NAN; }
    long double ExcessKurtosis() const override { return NAN; }

    /**
     * @fn CircularMean
     * @return Circular mean of random variable
     */
    virtual long double CircularMean() const = 0;

    /**
     * @fn CircularVariance
     * @return Circular variance of random variable
     */
    virtual long double CircularVariance() const = 0;
};

#endif // CIRCULARDISTRIBUTION_H

#ifndef CIRCULARDISTRIBUTION_H
#define CIRCULARDISTRIBUTION_H

#include "../ContinuousDistribution.h"

/**
 * @brief The CircularDistribution class <BR>
 * Abstract class for all continuous circular distributions
 *
 * Note that all the moments are now useless, we implement circular moments instead.
 */
class RANDLIBSHARED_EXPORT CircularDistribution : public ContinuousDistribution
{
protected:
    double loc{};

    CircularDistribution(double location = 0);
    virtual ~CircularDistribution() {}

public:
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return loc - M_PI; }
    double MaxValue() const override { return loc + M_PI; }

    void SetLocation(double location);
    inline double GetLocation() const { return loc; }

    double Mean() const override { return NAN; }
    double Variance() const override { return NAN; }
    double Skewness() const override { return NAN; }
    double ExcessKurtosis() const override { return NAN; }

    /**
     * @fn CircularMean
     * @return Circular mean of random variable
     */
    virtual double CircularMean() const = 0;

    /**
     * @fn CircularVariance
     * @return Circular variance of random variable
     */
    virtual double CircularVariance() const = 0;
};

#endif // CIRCULARDISTRIBUTION_H

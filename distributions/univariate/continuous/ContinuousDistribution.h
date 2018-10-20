#ifndef CONTINUOUS_DISTRIBUTION_H
#define CONTINUOUS_DISTRIBUTION_H

#include "../UnivariateDistribution.h"
#include "../ExponentialFamily.h"

/**
 * @brief The ContinuousDistribution class <BR>
 * Abstract class for all continuous distributions
 */
template < typename RealType >
class RANDLIBSHARED_EXPORT ContinuousDistribution : public virtual UnivariateDistribution<RealType>
{
    static_assert(std::is_floating_point_v<RealType>, "Continuous distribution supports only floating-point types");

protected:
    ContinuousDistribution() {}
    virtual ~ContinuousDistribution() {}

public:
    /**
     * @fn f
     * @param x
     * @return probability density function
     */
    virtual double f(const RealType & x) const { return std::exp(this->logf(x)); }

    /**
     * @fn logf
     * @param x
     * @return logarithm of probability density function
     */
    virtual double logf(const RealType & x) const = 0;

    /**
     * @fn ProbabilityDensityFunction
     * fill vector y with f(x)
     * @param x
     * @param y
     */
    void ProbabilityDensityFunction(const std::vector<RealType> &x, std::vector<double> &y) const;

    /**
     * @fn LogProbabilityDensityFunction
     * fill vector y with logf(x)
     * @param x
     * @param y
     */
    void LogProbabilityDensityFunction(const std::vector<RealType> &x, std::vector<double> &y) const;

    RealType Mode() const override;

protected:
    RealType quantileImpl(double p, RealType initValue) const override;
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p, RealType initValue) const override;
    RealType quantileImpl1m(double p) const override;
    long double ExpectedValue(const std::function<double (RealType)> &funPtr, RealType minPoint, RealType maxPoint) const override;

public:
    double Hazard(const RealType &x) const override;
    double LikelihoodFunction(const std::vector<RealType> &sample) const override;
    double LogLikelihoodFunction(const std::vector<RealType> &sample) const override;

    /**
     * @fn KolmogorovSmirnovTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha level of test
     * @return true if sample is from this distribution according to asymptotic KS-test, false otherwise
     */
    bool KolmogorovSmirnovTest(const std::vector<RealType> &orderStatistic, double alpha) const;
};


template < typename RealType, typename P >
class RANDLIBSHARED_EXPORT ContinuousExponentialFamily : public ExponentialFamily<RealType, P>, public ContinuousDistribution<RealType>
{
public:
    virtual double logf(const RealType & x) const {
        return this->LogProbabilityMeasure(x);
    }
};

#endif // CONTINUOUS_DISTRIBUTION_H

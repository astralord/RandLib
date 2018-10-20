#ifndef DISCRETE_DISTRIBUTION_H
#define DISCRETE_DISTRIBUTION_H

#include "../UnivariateDistribution.h"
#include "../ExponentialFamily.h"

/**
 *@brief The DiscreteDistribution class <BR>
 * Abstract class for all discrete distributions
 */
template < typename IntType >
class RANDLIBSHARED_EXPORT DiscreteDistribution : public virtual UnivariateDistribution<IntType>
{
    static_assert(std::is_integral_v<IntType> && std::is_signed_v<IntType>, "Discrete distribution supports only signed integral types");

protected:
    DiscreteDistribution() {}
    virtual ~DiscreteDistribution() {}

public:
    /**
     * @fn P
     * @param k
     * @return probability to get k
     */
    virtual double P(const IntType & k) const { return std::exp(this->logP(k)); }

    /**
     * @fn logP
     * @param x
     * @return logarithm of probability to get x
     */
    virtual double logP(const IntType & x) const = 0;

    /**
     * @fn ProbabilityMassFunction
     * fill vector y with P(x)
     * @param x
     * @param y
     */
    void ProbabilityMassFunction(const std::vector<IntType> &x, std::vector<double> &y) const;

    /**
     * @fn LogProbabilityMassFunction
     * fill vector y with logP(x)
     * @param x
     * @param y
     */
    void LogProbabilityMassFunction(const std::vector<IntType> &x, std::vector<double> &y) const;

    IntType Mode() const override;

protected:
    IntType quantileImpl(double p, IntType initValue) const override;
    IntType quantileImpl(double p) const override;
    IntType quantileImpl1m(double p, IntType initValue) const override;
    IntType quantileImpl1m(double p) const override;
    long double ExpectedValue(const std::function<double (IntType)> &funPtr, IntType minPoint, IntType maxPoint) const override;

public:
    /**
     * @fn Hazard
     * @param x
     * @return hazard function
     */
    double Hazard(const IntType &x) const override;

    /**
     * @fn LikelihoodFunction
     * @param sample
     * @return likelihood function of the distribution for given sample
     */
    double LikelihoodFunction(const std::vector<IntType> &sample) const override;

    /**
     * @fn LogLikelihoodFunction
     * @param sample
     * @return log-likelihood function of the distribution for given sample
     */
    double LogLikelihoodFunction(const std::vector<IntType> &sample) const override;

    /**
     * @fn PearsonChiSquaredTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha significance level of the test
     * @param lowerBoundary setting of the left most interval (-infinity, lowerBoundary]
     * @param upperBoundary setting of the right most interval [upperBoundary, infinity)
     * @param numberOfEstimatedParameters zero by default
     * @return true if sample is from this distribution according to Pearson's chi-squared test, false otherwise
     *
     * We wrote chi-squared test only for discrete distribution, as this is much more complicated for continuous one.
     * The problem is ambiguity of grouping sample into intervals. In the case when parameters are estimated by
     * maximum-likelihood estimator, using original observations, statistics might not follow asymptotic chi-square
     * distribution and that leads to serious underestimate of the error of the first kind.
     * For more details look: "The use of MLE in chi-square tests for goodness of fit" by Herman Chernoff and E.L. Lehmann
     */
    bool PearsonChiSquaredTest(const std::vector<IntType> &orderStatistic, double alpha, int lowerBoundary, int upperBoundary, size_t numberOfEstimatedParameters = 0) const;

    /**
     * @fn PearsonChiSquaredTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha significance level of the test
     * @param numberOfEstimatedParameters zero by default
     * @return true if sample is from this distribution according to Pearson's chi-squared test, false otherwise
     * In this function user won't set upper and lower intervals for tails.
     * However it might be useful to group rare events for chi-squared test to give better results
     */
    bool PearsonChiSquaredTest(const std::vector<IntType> &orderStatistic, double alpha, size_t numberOfEstimatedParameters = 0) const;
};



template < typename IntType, typename P >
class RANDLIBSHARED_EXPORT DiscreteExponentialFamily : public ExponentialFamily<IntType, P>, public DiscreteDistribution<IntType>
{
public:
    virtual double logP(const IntType & x) const {
        return this->LogProbabilityMeasure(x);
    }
};

#endif // DISCRETE_DISTRIBUTION_H

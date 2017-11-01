#ifndef DISCRETE_DISTRIBUTION_H
#define DISCRETE_DISTRIBUTION_H

#include "../UnivariateDistribution.h"

/**
 *@brief The DiscreteDistribution class <BR>
 * Abstract class for all discrete distributions
 */
class RANDLIBSHARED_EXPORT DiscreteDistribution : public virtual UnivariateDistribution<int>
{
protected:
    DiscreteDistribution() {}
    virtual ~DiscreteDistribution() {}

public:
    /**
     * @fn P
     * @param k
     * @return probability to get k
     */
    virtual double P(const int & k) const = 0;

    /**
     * @fn logP
     * @param x
     * @return logarithm of probability to get x
     */
    virtual double logP(const int & x) const = 0;

    /**
     * @fn ProbabilityMassFunction
     * fill vector y with P(x)
     * @param x
     * @param y
     */
    void ProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const;

    /**
     * @brief LogProbabilityMassFunction
     * fill vector y with logP(x)
     * @param x
     * @param y
     */
    void LogProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const;

    int Mode() const override;

private:
    int quantileImpl(double p) const override;
    int quantileImpl1m(double p) const override;
    double ExpectedValue(const std::function<double (double)> &funPtr, int minPoint, int maxPoint) const override;

public:
    double Hazard(double x) const override;

    /**
     * @brief LikelihoodFunction
     * @param sample
     * @return likelihood function of the distribution for given sample
     */
    double LikelihoodFunction(const std::vector<int> &sample) const override;

    /**
     * @brief LogLikelihoodFunction
     * @param sample
     * @return log-likelihood function of the distribution for given sample
     */
    double LogLikelihoodFunction(const std::vector<int> &sample) const override;

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
    bool PearsonChiSquaredTest(const std::vector<int> &orderStatistic, double alpha, int lowerBoundary, int upperBoundary, size_t numberOfEstimatedParameters = 0) const;

    /**
     * @brief PearsonChiSquaredTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha significance level of the test
     * @param numberOfEstimatedParameters zero by default
     * @return true if sample is from this distribution according to Pearson's chi-squared test, false otherwise
     * In this function user won't set upper and lower intervals for tails.
     * However it might be useful to group rare events for chi-squared test to give better results
     */
    bool PearsonChiSquaredTest(const std::vector<int> &orderStatistic, double alpha, size_t numberOfEstimatedParameters = 0) const;
};

#endif // DISCRETE_DISTRIBUTION_H

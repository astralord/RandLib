#ifndef DISCRETE_DISTRIBUTION_H
#define DISCRETE_DISTRIBUTION_H

#include "../UnivariateProbabilityDistribution.h"

/**
 *@brief The DiscreteDistribution class <BR>
 * Abstract class for all discrete distributions
 */
class RANDLIBSHARED_EXPORT DiscreteDistribution : public virtual UnivariateProbabilityDistribution<int>
{
public:
    DiscreteDistribution() {}
    virtual ~DiscreteDistribution() {}

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

    int Mode() const override;

private:
    int quantileImpl(double p) const override;
    int quantileImpl1m(double p) const override;
    double ExpectedValue(const std::function<double (double)> &funPtr, int minPoint, int maxPoint) const override;

public:
    double Hazard(double x) const override;

    double Likelihood(const std::vector<int> &sample) const override;
    double LogLikelihood(const std::vector<int> &sample) const override;

    /**
     * @fn PearsonChiSquaredTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha level of test
     * @param numberOfEstimatedParameters zero by default
     * @param lowerBoundary setting of the left most interval (-infinity, lowerBoundary]
     * @param upperBoundary setting of the right most interval [upperBoundary, infinity)
     * @param numberOfEstimatedParameters
     * @return true if sample is from this distribution according to Pearson's chi-squared test, false otherwise
     */
    bool PearsonChiSquaredTest(const std::vector<int> &orderStatistic, double alpha, int lowerBoundary, int upperBoundary, size_t numberOfEstimatedParameters = 0) const;
    bool PearsonChiSquaredTest(const std::vector<int> &orderStatistic, double alpha, size_t numberOfEstimatedParameters = 0) const;
};

#endif // DISCRETE_DISTRIBUTION_H

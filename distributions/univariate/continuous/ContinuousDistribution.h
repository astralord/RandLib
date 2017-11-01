#ifndef CONTINUOUS_DISTRIBUTION_H
#define CONTINUOUS_DISTRIBUTION_H

#include "../UnivariateDistribution.h"

/**
 * @brief The ContinuousDistribution class <BR>
 * Abstract class for all continuous distributions
 */
class RANDLIBSHARED_EXPORT ContinuousDistribution : public virtual UnivariateDistribution<double>
{
protected:
    ContinuousDistribution() {}
    virtual ~ContinuousDistribution() {}

public:
    /**
     * @fn f
     * @param x
     * @return probability density function
     */
    virtual double f(const double & x) const = 0;

    /**
     * @fn logf
     * @param x
     * @return logarithm of probability density function
     */
    virtual double logf(const double & x) const = 0;

    /**
     * @fn ProbabilityDensityFunction
     * fill vector y by f(x)
     * @param x
     * @param y
     */
    void ProbabilityDensityFunction(const std::vector<double> &x, std::vector<double> &y) const;

    /**
     * @brief LogProbabilityDensityFunction
     * fill vector y by logf(x)
     * @param x
     * @param y
     */
    void LogProbabilityDensityFunction(const std::vector<double> &x, std::vector<double> &y) const;

    double Mode() const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
    double ExpectedValue(const std::function<double (double)> &funPtr, double minPoint, double maxPoint) const override;

public:
    double Hazard(double x) const override;
    double LikelihoodFunction(const std::vector<double> &sample) const override;
    double LogLikelihoodFunction(const std::vector<double> &sample) const override;

    /**
     * @fn KolmogorovSmirnovTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha level of test
     * @return true if sample is from this distribution according to asymptotic KS-test, false otherwise
     */
    bool KolmogorovSmirnovTest(const std::vector<double> &orderStatistic, double alpha) const;
};

#endif // CONTINUOUS_DISTRIBUTION_H

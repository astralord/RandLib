#ifndef CONTINUOUS_DISTRIBUTION_H
#define CONTINUOUS_DISTRIBUTION_H

#include "../UnivariateProbabilityDistribution.h"

/**
 * @brief The ContinuousDistribution class
 */
class RANDLIBSHARED_EXPORT ContinuousDistribution : public virtual UnivariateProbabilityDistribution<double>
{
public:
    ContinuousDistribution() {}
    virtual ~ContinuousDistribution() {}

    /**
     * @brief f
     * probability density function
     * @param x
     * @return
     */
    virtual double f(double x) const = 0;

    /**
     * @brief ProbabilityDensityFunction
     * fill vector y by f(x)
     * @param x
     * @param y
     */
    void ProbabilityDensityFunction(const std::vector<double> &x, std::vector<double> &y) const;

    double Mode() const override;

protected:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;
    double ExpectedValue(const std::function<double (double)> &funPtr, double minPoint, double maxPoint) const override;

public:
    double Hazard(double x) const override;
    double Likelihood(const std::vector<double> &sample) const override;
    double LogLikelihood(const std::vector<double> &sample) const override;

    /**
     * @brief KolmogorovSmirnovTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha level of test
     * @return true if sample is from this distribution according to asymptotic KS-test, false otherwise
     */
    bool KolmogorovSmirnovTest(const std::vector<double> &orderStatistic, double alpha) const;
};

#endif // CONTINUOUS_DISTRIBUTION_H

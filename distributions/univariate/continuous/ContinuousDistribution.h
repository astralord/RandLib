#ifndef CONTINUOUS_DISTRIBUTION_H
#define CONTINUOUS_DISTRIBUTION_H

#include "../UnivariateProbabilityDistribution.h"

/**
 * @brief The ContinuousDistribution class
 */
class RANDLIBSHARED_EXPORT ContinuousDistribution : public UnivariateProbabilityDistribution
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
     * @brief pdf
     * fill vector y by f(x)
     * @param x
     * @param y
     */
    void pdf(const std::vector<double> &x, std::vector<double> &y) const;

    double Quantile(double p) const override;
    double Hazard(double x) const override;

    double Median() const override;
    double Mode() const override;

protected:
    double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const override;

public:
    double likelihood(const std::vector<double> &sample) const;
    double loglikelihood(const std::vector<double> &sample) const;
};

#endif // CONTINUOUS_DISTRIBUTION_H

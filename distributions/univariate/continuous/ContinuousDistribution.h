#ifndef CONTINUOUS_DISTRIBUTION_H
#define CONTINUOUS_DISTRIBUTION_H

#include "../UnivariateProbabilityDistribution.h"

template <typename T>
class SumContinuousRand;

/**
 * @brief The ContinuousDistribution class
 */
class RANDLIBSHARED_EXPORT ContinuousDistribution : public UnivariateProbabilityDistribution<double>
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

    double Quantile(double p) const override;
    double Hazard(double x) const override;
    double Mode() const override;

protected:
    double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const override;

public:
    double Likelihood(const std::vector<double> &sample) const override;
    double LogLikelihood(const std::vector<double> &sample) const override;

    friend const SumContinuousRand<double> operator+(const ContinuousDistribution& left, const ContinuousDistribution& right);
};

#endif // CONTINUOUS_DISTRIBUTION_H

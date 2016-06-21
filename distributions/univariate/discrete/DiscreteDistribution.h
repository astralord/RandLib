#ifndef DISCRETE_DISTRIBUTION_H
#define DISCRETE_DISTRIBUTION_H

#include "../UnivariateProbabilityDistribution.h"

/**
 *@brief The DiscreteDistribution class
 */
class RANDLIBSHARED_EXPORT DiscreteDistribution : public virtual UnivariateProbabilityDistribution<int>
{
public:
    DiscreteDistribution() {}
    virtual ~DiscreteDistribution() {}

    /**
     * @brief P
     * probability to get x
     * @param x
     * @return
     */
    virtual double P(int x) const = 0;

    /**
     * @brief ProbabilityMassFunction
     * fill vector y by P(x)
     * @param x
     * @param y
     */
    void ProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const;

    double Quantile(double probability) const override;

    double Hazard(double x) const override;

public:
    double ExpectedValue(const std::function<double (double)> &funPtr, int minPoint, int maxPoint) const override;
    double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const override;

    double Likelihood(const std::vector<int> &sample) const override;
    double LogLikelihood(const std::vector<int> &sample) const override;
};

#endif // DISCRETE_DISTRIBUTION_H

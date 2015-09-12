#ifndef CONTINUOUSRAND_H
#define CONTINUOUSRAND_H

#include "../RandomVariable.h"

/**
 * @brief The ContinuousRand class
 */
class RANDLIBSHARED_EXPORT ContinuousRand : public RandomVariable
{
public:
    ContinuousRand() : RandomVariable() {}
    virtual ~ContinuousRand() {}

    /**
     * @brief f
     * probability density function
     * @param x
     * @return
     */
    virtual double f(double x) const = 0;

    void pdf(const QVector<double> &x, QVector<double> &y) const;

    double Quantile(double p) const override;

    double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const override;
    double Median() const override;
    double Mode() const override;

    double likelihood(const QVector<double> &sample) const;
    double loglikelihood(const QVector<double> &sample) const;
};

#endif // CONTINUOUSRAND_H

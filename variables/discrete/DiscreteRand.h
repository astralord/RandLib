#ifndef DISCRETERAND_H
#define DISCRETERAND_H

#include "../RandomVariable.h"

/**
 *@brief The DiscreteRand class
 */
class RANDLIBSHARED_EXPORT DiscreteRand : public RandomVariable
{
public:
    DiscreteRand() : RandomVariable() {}
    virtual ~DiscreteRand() {}

    /**
     * @brief P
     * probability to get x
     * @param x
     * @return
     */
    virtual double P(int x) const = 0;

    void pmf(const QVector<int> &x, QVector<double> &y) const;

    double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const override;

    double likelihood(const QVector<int> &sample) const;
    double loglikelihood(const QVector<int> &sample) const;
};

#endif // DISCRETERAND_H

#ifndef DISCRETERAND_H
#define DISCRETERAND_H

#include "../RandomVariable.h"

/**
 *@brief The DiscreteRand class
 */
template< typename T >
class RANDLIBSHARED_EXPORT DiscreteRand : public RandomVariable
{
public:
    DiscreteRand() : RandomVariable() {}
    virtual ~DiscreteRand() {}

    virtual double P(T x) const = 0;

    void pmf(const QVector<T> &x, QVector<double> &y) const;

    double likelihood(const QVector<T> &sample) const;
    double loglikelihood(const QVector<T> &sample) const;
};

#endif // DISCRETERAND_H

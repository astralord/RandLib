#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <QVector>
#include <QtMath>
#include <RandMath.h>
#include "randlib_global.h"

/**
 * @brief The RandomVariable class
 */
class RANDLIBSHARED_EXPORT RandomVariable {

    /**
     * @brief randValue
     * real value of random variable
     * (totally private)
     */
    unsigned long randValue;

protected:

    //TODO: add other algorithms
    unsigned long SHR3();

public:
    RandomVariable();
    void seedRand(unsigned long seed);

    virtual double cdf(double x) = 0;
    virtual double value() = 0;

    virtual double M() = 0;
    virtual double Var() = 0;

    void _cdf(const std::vector<double> &x, std::vector<double> &y);
    void sample(std::vector<double> &outputData);

    void _cdf(const QVector<double> &x, QVector<double> &y);
    void sample(QVector<double> &outputData);
};


/**
 * @brief The ContinuousRandom class
 */
class RANDLIBSHARED_EXPORT ContinuousRand : public RandomVariable
{
public:
    ContinuousRand() : RandomVariable() {}

    virtual double pdf(double x) = 0;

    void _pdf(const std::vector<double> &x, std::vector<double> &y);
    void _pdf(const QVector<double> &x, QVector<double> &y);

    double likelihood(const std::vector<double> &value);
    double loglikelihood(const std::vector<double> &value);
};

/**
 * @brief The DiscreteRandom class
 */

template <typename T>
class RANDLIBSHARED_EXPORT DiscreteRand : public RandomVariable
{
public:
    virtual double P(T x) = 0;
};

#endif // RANDOMVARIABLE_H

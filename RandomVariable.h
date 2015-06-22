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

    virtual double cdf(double x) const = 0;
    virtual double value() = 0;

    virtual double M() const = 0;
    virtual double Var() const = 0;

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

    virtual double pdf(double x) const = 0;

    void _pdf(const std::vector<double> &x, std::vector<double> &y) const;
    void _pdf(const QVector<double> &x, QVector<double> &y) const;

    double likelihood(const std::vector<double> &value) const;
    double loglikelihood(const std::vector<double> &value) const;
};



class RANDLIBSHARED_EXPORT DiscreteIntRand : public RandomVariable
{
public:
    virtual double P(int x) const = 0;

    void probs(const std::vector<int> &x, std::vector<double> &y) const;
    void probs(const QVector<int> &x, QVector<double> &y) const;
};

class RANDLIBSHARED_EXPORT DiscreteDoubleRand : public RandomVariable
{
public:
    virtual double P(double x) const = 0;

    void probs(const std::vector<double> &x, std::vector<double> &y) const;
    void probs(const QVector<double> &x, QVector<double> &y) const;
};

#endif // RANDOMVARIABLE_H

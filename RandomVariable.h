#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <QVector>
#include <QtMath>
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
    unsigned long randValue; // for SHR3

protected:
    /// Constants for optimization
    static constexpr double M_1_SQRTPI = 0.56418958354775628695;
    static constexpr double M_1_SQRT2PI = 0.39894228040143267794;

    static constexpr double MIN_POSITIVE = 1e-21;

    //TODO: add other algorithms
    unsigned long SHR3();

    /// useful functions
    // TODO: create math class for it
    double factorial(int n);
    double doubleFactorial(int n);
<<<<<<< HEAD
    double lowerIncGamma(double a, double x);
    double upperIncGamma(double a, double x);
=======
>>>>>>> 40068e8485c8670ce8d7f39d6c822f820d31fc77

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

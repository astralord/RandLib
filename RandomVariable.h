#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <vector>
#include <QVector>
#include <QtMath>

#include "randlib_global.h"

/**
 * @brief The RandomVariable class
 */
class RANDLIBSHARED_EXPORT RandomVariable
{
    /**
     * @brief randValue
     * real value of random variable
     * (totally private)
     */
    unsigned long randValue; // for SHR3
    //unsigned long long z, w; // for MWC

protected:
    /// Constants for optimization
    static constexpr double M_1_SQRTPI = 0.56418958354775628695;
    static constexpr double M_1_SQRT2PI = 0.39894228040143267794;

    static constexpr double MIN_POSITIVE = 1e-21;

    //TODO: add other algorithms
    unsigned long SHR3();
    //MWC doesn't work!
    unsigned long long MWC();

public:
    RandomVariable();
    void seedRand(unsigned long seed);

    virtual double pdf(double x) = 0;
    virtual double cdf(double x) = 0;
    virtual double value() = 0;

    virtual double M() = 0;
    virtual double Var() = 0;

    void _pdf(const std::vector<double> &x, std::vector<double> &y);
    void _cdf(const std::vector<double> &x, std::vector<double> &y);
    void sample(std::vector<double> &outputData);
    void _pdf(const QVector<double> &x, QVector<double> &y);
    void _cdf(const QVector<double> &x, QVector<double> &y);
    void sample(QVector<double> &outputData);

    double likelihood(const std::vector<double> &value);
    double loglikelihood(const std::vector<double> &value);

};

#endif // RANDOMVARIABLE_H

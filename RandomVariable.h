#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <QVector>
#include <QtMath>
#include <RandMath.h>
#include "randlib_global.h"

//TODO: implement
// Beta ++-
// Gamma -++
// Stable +--
// Fisher Snedecor ++-
// Binomial ---
// Multinomial ---
// Weibull ---
// Triangular ---
// Logistic ---
// Rayleigh ---
// ALL: http://www.mathworks.com/help/stats/makedist.html

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
    virtual ~RandomVariable() {}

    void seedRand(unsigned long seed);

    virtual double F(double x) const = 0;

    virtual double value() = 0;
    void sample(std::vector<double> &outputData);
    void sample(QVector<double> &outputData);

    virtual double M() const = 0;
    virtual double Var() const = 0;

};


/**
 * @brief The ContinuousRand class
 */
class RANDLIBSHARED_EXPORT ContinuousRand : public RandomVariable
{
public:
    ContinuousRand() : RandomVariable() {}
    virtual ~ContinuousRand() {}

    virtual double f(double x) const = 0;

    void pdf(const std::vector<double> &x, std::vector<double> &y) const;
    void pdf(const QVector<double> &x, QVector<double> &y) const;

    void cdf(const std::vector<double> &x, std::vector<double> &y);
    void cdf(const QVector<double> &x, QVector<double> &y);

    double likelihood(const std::vector<double> &sample) const;
    double loglikelihood(const std::vector<double> &sample) const;
};


//TODO: make value integer!
class RANDLIBSHARED_EXPORT DiscreteIntRand : public RandomVariable
{
public:
    DiscreteIntRand() : RandomVariable() {}
    virtual ~DiscreteIntRand() {}

    virtual double P(int x) const = 0;

    void pmf(const std::vector<int> &x, std::vector<double> &y) const;
    void pmf(const QVector<int> &x, QVector<double> &y) const;

    void cdf(const std::vector<int> &x, std::vector<double> &y);
    void cdf(const QVector<int> &x, QVector<double> &y);

    double likelihood(const std::vector<int> &sample) const;
    double loglikelihood(const std::vector<int> &sample) const;
};


class RANDLIBSHARED_EXPORT DiscreteDoubleRand : public RandomVariable
{
public:
    DiscreteDoubleRand() : RandomVariable() {}
    virtual ~DiscreteDoubleRand() {}

    virtual double P(double x) const = 0;

    void pmf(const std::vector<double> &x, std::vector<double> &y) const;
    void pmf(const QVector<double> &x, QVector<double> &y) const;

    void cdf(const std::vector<double> &x, std::vector<double> &y);
    void cdf(const QVector<double> &x, QVector<double> &y);

    double likelihood(const std::vector<double> &sample) const;
    double loglikelihood(const std::vector<double> &sample) const;
};

#endif // RANDOMVARIABLE_H

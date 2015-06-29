#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <QVector>
#include <QtMath>
#include "RandMath.h"
#include "randlib_global.h"

//TODO: Distribution (value, pdf, cdf)
// Beta (++-)
// Stable (+--)
// Fisher Snedecor (++-)
// StudentTRand (++-)

// Binomial (---)
// Multinomial (---)
// Weibull (---)
// Triangular (---)
// Logistic (---)
// Rayleigh (---)
// ALL: http://www.mathworks.com/help/stats/makedist.html

/**
 * @brief The RandomVariable class
 */
template< typename T >
class RANDLIBSHARED_EXPORT RandomVariable {

    /// Variables for pseudo generator
    static unsigned long startPoint;
    unsigned long X, Y, Z, W;
    bool C;

    /// Variables for quasi generator
    unsigned N, CODE, Q;
protected:
    /**
     * @brief fastKISS
     * Random variable generator Keep It Simply Stupid
     * @return random variable uniformly distributed from 0 to 2^32-1
     */
    unsigned long fastKISS();

    unsigned long quasiGen();

public:
    RandomVariable();
    virtual ~RandomVariable() {}

    virtual T value() = 0;

    void sample(QVector<T> &outputData);

    /**
     * @brief F
     * @param x
     * @return P(X < x)
     */
    virtual double F(double x) const = 0;

    void cdf(const QVector<double> &x, QVector<double> &y);

    /**
     * @brief M
     * @return Mathematical expectation
     */
    virtual double M() const = 0;
    /**
     * @brief Var
     * @return Variance
     */
    virtual double Var() const = 0;
};


/**
 * @brief The ContinuousRand class
 */
class RANDLIBSHARED_EXPORT ContinuousRand : public RandomVariable<double>
{
public:
    ContinuousRand() : RandomVariable() {}
    virtual ~ContinuousRand() {}

    virtual double f(double x) const = 0;

    void pdf(const QVector<double> &x, QVector<double> &y) const;

    double likelihood(const QVector<double> &sample) const;
    double loglikelihood(const QVector<double> &sample) const;
};

/**
 *@brief The DiscreteRand class
 */
template< typename T >
class RANDLIBSHARED_EXPORT DiscreteRand : public RandomVariable<T>
{
public:
    DiscreteRand() : RandomVariable<T>() {}
    virtual ~DiscreteRand() {}

    virtual double P(T x) const = 0;

    void pmf(const QVector<T> &x, QVector<double> &y) const;

    double likelihood(const QVector<T> &sample) const;
    double loglikelihood(const QVector<T> &sample) const;
};

#endif // RANDOMVARIABLE_H

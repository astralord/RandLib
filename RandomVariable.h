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
class RANDLIBSHARED_EXPORT RandomVariable {

    static unsigned long startPoint;
    unsigned long X, Y, Z, W;
    bool C;

protected:
    /**
     * @brief fastKISS
     * Random variable generator Keep It Simply Stupid
     * @return random variable uniformly distributed from 0 to 2^32-1
     */
    unsigned long fastKISS();

public:
    RandomVariable();
    virtual ~RandomVariable() {}
    /**
     * @brief F
     * @param x
     * @return P(X < x)
     */
    virtual double F(double x) const = 0;
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
 * @brief The IntRand class
 * Distribution over integer values
 * (beware: values cannot be greater than 2^32)
 */
class RANDLIBSHARED_EXPORT IntRand : public RandomVariable
{
public:
    IntRand() : RandomVariable() {}
    virtual ~IntRand() {}

    virtual int value() = 0;

    void sample(std::vector<int> &outputData);
    void sample(QVector<int> &outputData);
};

/**
 * @brief The DoubleRand class
 */
class RANDLIBSHARED_EXPORT DoubleRand : public RandomVariable
{
public:
    DoubleRand() : RandomVariable() {}
    virtual ~DoubleRand() {}

    virtual double value() = 0;

    void sample(std::vector<double> &outputData);
    void sample(QVector<double> &outputData);
};

/**
 * @brief The ContinuousRand class
 */
class RANDLIBSHARED_EXPORT ContinuousRand : public DoubleRand
{
public:
    ContinuousRand() : DoubleRand() {}
    virtual ~ContinuousRand() {}

    virtual double f(double x) const = 0;

    void pdf(const std::vector<double> &x, std::vector<double> &y) const;
    void pdf(const QVector<double> &x, QVector<double> &y) const;

    void cdf(const std::vector<double> &x, std::vector<double> &y);
    void cdf(const QVector<double> &x, QVector<double> &y);

    double likelihood(const std::vector<double> &sample) const;
    double loglikelihood(const std::vector<double> &sample) const;
};


/**
 * @brief The DiscreteIntRand class
 */
class RANDLIBSHARED_EXPORT DiscreteIntRand : public IntRand
{
public:
    DiscreteIntRand() : IntRand() {}
    virtual ~DiscreteIntRand() {}

    virtual double P(int x) const = 0;

    void pmf(const std::vector<int> &x, std::vector<double> &y) const;
    void pmf(const QVector<int> &x, QVector<double> &y) const;

    void cdf(const std::vector<int> &x, std::vector<double> &y) const;
    void cdf(const QVector<int> &x, QVector<double> &y) const;

    double likelihood(const std::vector<int> &sample) const;
    double loglikelihood(const std::vector<int> &sample) const;
};


/**
 * @brief The DiscreteDoubleRand class
 */
class RANDLIBSHARED_EXPORT DiscreteDoubleRand : public DoubleRand
{
public:
    DiscreteDoubleRand() : DoubleRand() {}
    virtual ~DiscreteDoubleRand() {}

    virtual double P(double x) const = 0;

    void pmf(const std::vector<double> &x, std::vector<double> &y) const;
    void pmf(const QVector<double> &x, QVector<double> &y) const;

    void cdf(const std::vector<double> &x, std::vector<double> &y) const;
    void cdf(const QVector<double> &x, QVector<double> &y) const;

    double likelihood(const std::vector<double> &sample) const;
    double loglikelihood(const std::vector<double> &sample) const;
};

#endif // RANDOMVARIABLE_H

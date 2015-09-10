#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <cmath>
#include <string>

#include "math/RandMath.h"
#include "randlib_global.h"

#include <QDebug>
#include <QVector>

/**
 * @brief The RandomVariable class
 */
class RANDLIBSHARED_EXPORT RandomVariable
{

protected:

    std::string toStringWithPrecision(const double a_value, const int n = 6);

public:
    RandomVariable();
    virtual ~RandomVariable() {}

    /**
     * @brief name
     * @return name of distribution, for instance "Normal(0, 1)"
     */
    virtual std::string name() = 0;

    /**
     * @brief variate()
     * @return random variable
     */
    virtual double variate() const = 0;

    /**
     * @brief sample
     * @param outputData
     */
    virtual void sample(QVector<double> &outputData);

    /**
     * @brief F
     * @param x
     * @return P(X < x)
     */
    virtual double F(double x) const = 0;

    /**
     * @brief cdf
     * @param x input vector
     * @param y output vector: y = P(X < x)
     */
    void cdf(const QVector<double> &x, QVector<double> &y);

    /**
     * @brief M
     * @return Mathematical expectation
     */
    virtual double E() const = 0;

    /**
     * @brief Var
     * @return Variance
     */
    virtual double Var() const = 0;

    /**
     * @brief CF
     * @param x
     * @return Characteristic function (inverse Fourier transform of probability function)
     */
    virtual std::complex<double> CF(double t) const // = 0
    {
        return std::complex<double>(t);
    }

    /**
     * @brief cf
     * @param x input vector
     * @param y output vector: y = CF(x)
     */
    void cf(const QVector<double> &t, QVector<std::complex<double>> &y);

    /**
     * @brief quantile
     * @param p
     * @return such x that F(x) = p
     */
    virtual double quantile(double p) const {
        if (p < 0 || p > 1)
            return NAN;
        double root = 0;
        RandMath::findRoot([this, p] (double x)
        {
            return F(x) - p;
        },
        0, 1, root);
        return root;
    }

    /**
     * @brief Median
     * @return such x that F(x) = 0.5
     */
    virtual double Median() const {
        return quantile(0.5);
    }

    /**
     * @brief Mode
     * @return the most probable value
     */
    virtual double Mode() const { return 0; }

    /**
     * @brief Skewness
     * @return E[((X - mu) / sigma) ^ 3]
     * where mu is central moment and sigma is standard deviation
     */
    virtual double Skewness() const { return 0; }

    /**
     * @brief ExcessKurtosis
     * @return mu_4 / sigma ^ 4 - 3
     * (fourth moment around the mean divided by the square of the variance of the probability distribution minus 3)
     */
    virtual double ExcessKurtosis() const { return 0; }
};

#endif // RANDOMVARIABLE_H

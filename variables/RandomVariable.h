#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <cmath>
#include <string>

#include "math/RandMath.h"
#include "randlib_global.h"

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
    virtual void sample(QVector<double> &outputData) const;

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
    void cdf(const QVector<double> &x, QVector<double> &y) const;

    /**
     * @brief Mean
     * @return Mathematical expectation
     */
    virtual double Mean() const = 0;

    /**
     * @brief Variance
     * @return Variance of random variable
     */
    virtual double Variance() const = 0;
    
    /**
     * @brief quantile
     * @param p
     * @return such x that F(x) = p
     */
    virtual double Quantile(double p) const = 0;
    
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
    void cf(const QVector<double> &t, QVector<std::complex<double>> &y) const;

    /**
     * @brief Hazard
     * return hazard function: pdf (or pmf) / (1 - cdf)
     * @param x input parameter
     */
    virtual double Hazard(double x) const = 0;
    
    /**
     * @brief ExpectedValue
     * @param funPtr function which expected value should be returned
     * @param startPoint argument in which vicinity value of funPtr definitely wouldn't be zero
     * @return E[funPtr(x)]
     */
    virtual double ExpectedValue(const std::function<double (double)> &funPtr, double startPoint) const = 0;

    /**
     * @brief Median
     * @return such x that F(x) = 0.5
     */
    virtual double Median() const;

    /**
     * @brief Mode
     * @return the most probable value
     */
    virtual double Mode() const = 0;

    /**
     * @brief Skewness
     * @return E[((X - mu) / sigma) ^ 3]
     * where mu is central moment and sigma is standard deviation
     */
    virtual double Skewness() const;

    /**
     * @brief ExcessKurtosis
     * @return E[((X - mu) / sigma) ^ 4]  - 3
     * (fourth moment around the mean divided by the square of the variance of the probability distribution minus 3)
     */
    virtual double ExcessKurtosis() const;

    /**
     * @brief Kurtosis
     * @return unbiased kurtosis = mu_4 / sigma ^ 4
     */
    virtual double Kurtosis() const;
};

#endif // RANDOMVARIABLE_H

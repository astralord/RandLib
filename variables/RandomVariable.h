#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <QVector>
#include <QtMath>
#include "math/RandMath.h"
#include "randlib_global.h"

#include <QDebug>

/**
 * @brief The RandomVariable class
 */
template< typename T >
class RANDLIBSHARED_EXPORT RandomVariable
{

public:
    RandomVariable();
    virtual ~RandomVariable() {}

    /**
     * @brief variate()
     * @return random variable
     */
    virtual T variate() = 0;

    /**
     * @brief sample
     * @param outputData
     */
    void sample(QVector<T> &outputData);

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
};

#endif // RANDOMVARIABLE_H

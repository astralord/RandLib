#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include <QVector>
#include <QtMath>
#include <string>
#include "math/RandMath.h"
#include "randlib_global.h"

#include <QDebug>

/**
 * @brief The RandomVariable class
 */
class RANDLIBSHARED_EXPORT RandomVariable
{

protected:
    std::string nameStr;

    std::string toStringWithPrecision(const double a_value, const int n = 6);

public:
    RandomVariable();
    virtual ~RandomVariable() {}

    /**
     * @brief name
     * @return
     */
    std::string name() const { return nameStr; }

    /**
     * @brief setName
     */
    virtual void setName() = 0;

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
};

#endif // RANDOMVARIABLE_H

#ifndef STOCHASTICPROCESS_H
#define STOCHASTICPROCESS_H

#include "randlib_global.h"
#include <QVector>

/**
 * @brief The StochasticProcess class
 */
class RANDLIBSHARED_EXPORT StochasticProcess
{
public:
    StochasticProcess();

    /**
     * @brief Mean
     * @return Mathematical expectation
     */
    virtual void Mean(const QVector<double> &time, QVector<double> &output) const = 0;
    /**
     * @brief Variance
     * @return Variance = E[X(t)^2] - E[X(t)]^2
     */
    virtual void Variance(const QVector<double> &time, QVector<double> &output) const = 0;
};

#endif // STOCHASTICPROCESS_H

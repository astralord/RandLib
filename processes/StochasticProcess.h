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
     * @brief M
     * @return Mathematical expectation
     */
    virtual void E(const QVector<double> &time, QVector<double> &output) const = 0;
    /**
     * @brief Var
     * @return Variance
     */
    virtual void Var(const QVector<double> &time, QVector<double> &output) const = 0;
};

#endif // STOCHASTICPROCESS_H

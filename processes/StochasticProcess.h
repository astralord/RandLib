#ifndef STOCHASTICPROCESS_H
#define STOCHASTICPROCESS_H

#include "randlib_global.h"
#include <QVector>

/**
 * @brief The StochasticProcess class
 */
class RANDLIBSHARED_EXPORT StochasticProcess
{
protected:
    double currentTime, currentValue, dt;

    virtual double nextImpl() = 0;
    virtual double nextImpl(double deltaT) = 0;
    virtual double MeanImpl(double t) const = 0;
    virtual double VarianceImpl(double t) const = 0;
    virtual double QuantileImpl(double t, double p) const = 0;

public:
    explicit StochasticProcess(double deltaT = 1.0);

    inline double getCurrentTime() const { return currentTime; }
    inline double getCurrentValue() const { return currentValue; }

    /**
     * @brief next
     * @return
     */
    double next();

    /**
     * @brief next
     * @param deltaT
     * @return
     */
    double next(double deltaT);

    /**
     * @brief Mean
     * @param t
     * @return Conditional mathematical expectation E[X(t)|F_s]
     */
    double Mean(double t) const;

    /**
     * @brief Variance
     * @param t
     * @return Conditional variance Var(X(t)|F_s)
     */
    double Variance(double t) const;

    /**
     * @brief Quantile
     * @param t
     * @param p
     * @return
     */
    double Quantile(double t, double p) const;
};

#endif // STOCHASTICPROCESS_H

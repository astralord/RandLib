#ifndef STOCHASTICPROCESS_H
#define STOCHASTICPROCESS_H

#include "randlib_global.h"
#include "../math/RandMath.h"


/**
 * @brief The StochasticProcess class
 */
template <typename T>
class RANDLIBSHARED_EXPORT StochasticProcess
{
protected:
    double currentTime, dt;
    T currentValue;

public:
    explicit StochasticProcess(double deltaT = 1.0, T initialValue = 0);

    /**
     * @brief reset
     * @param initialValue
     * set time to 0 and current value to initialValue
     */
    void reset(T initialValue = 0.0);

    /**
     * @brief getCurrentTime
     * @return current time
     */
    inline double getCurrentTime() const { return currentTime; }

    /**
     * @brief getCurrentValue
     * @return current value
     */
    inline T getCurrentValue() const { return currentValue; }

    /**
     * @brief next
     * @return next value
     */
    T next();

    /**
     * @brief Mean
     * @param t
     * @return Conditional mathematical expectation E[X(t)|X(s)], where s is the current time
     */
    double Mean(double t) const;

    /**
     * @brief Variance
     * @param t
     * @return Conditional variance Var(X(t)|X(s)), where s is the current time
     */
    double Variance(double t) const;

private:
    virtual void nextImpl() = 0;
    virtual double MeanImpl(double t) const = 0;
    virtual double VarianceImpl(double t) const = 0;
};

#endif // STOCHASTICPROCESS_H

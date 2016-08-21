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
     * @brief reSet
     * @param initialValue
     * Set time to 0 and current value to initialValue
     */
    void reSet(T initialValue = 0.0);

    /**
     * @brief GetCurrentTime
     * @return current time
     */
    inline double GetCurrentTime() const { return currentTime; }

    /**
     * @brief GetCurrentValue
     * @return current value
     */
    inline T GetCurrentValue() const { return currentValue; }

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

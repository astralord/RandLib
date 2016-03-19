#ifndef UHLENBECKORNSTEINPROCESS_H
#define UHLENBECKORNSTEINPROCESS_H

#include "StochasticProcess.h"
#include "WienerProcess.h"
#include "../distributions//continuous/NormalRand.h"

class RANDLIBSHARED_EXPORT UhlenbeckOrnsteinProcess : public StochasticProcess
{
    double mu, sigma, theta, x0;
    WienerProcess W;

    mutable double lastRateTime, lastExpRateTime;
public:
    UhlenbeckOrnsteinProcess(double deltaT = 1.0, double mean = 0.0, double volatility = 1.0, double rate = 1.0, double initialValue = 0.0);

    void setParameters(double mean, double volatility, double rate, double initialValue);

    inline double getMean() { return mu; }
    inline double getVolatility() { return sigma; }
    inline double getRate() { return theta; }
    inline double getInitialValue() { return x0; }

    double next() const;
    double next(double deltaT) const;

    void Mean(const QVector<double> &time, QVector<double> &output) const;
    void Variance(const QVector<double> &time, QVector<double> &output) const;
};

#endif // UHLENBECKORNSTEINPROCESS_H

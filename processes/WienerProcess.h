#ifndef WIENERPROCESS_H
#define WIENERPROCESS_H

#include "StochasticProcess.h"
#include "../variables/continuous/NormalRand.h"

/**
 * @brief The WienerProcess class
 */
class RANDLIBSHARED_EXPORT WienerProcess : public StochasticProcess
{
    double mu, var, sigma;
    double muT, sigmaT;

    mutable double lastValue;
public:
    WienerProcess(double deltaT = 1.0, double mean = 0, double variance = 1);
    void setMean(double mean);
    void setVar(double variance);
    void setSigma(double volatility);

    inline double getMean() { return mu; }
    inline double getVar() { return var; }
    inline double getSigma() { return sigma; }

    double next() const;
    double next(double deltaT) const;

    void Mean(const QVector<double> &time, QVector<double> &output) const;
    void Variance(const QVector<double> &time, QVector<double> &output) const;
};

#endif // WIENERPROCESS_H

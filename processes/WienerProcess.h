#ifndef WIENERPROCESS_H
#define WIENERPROCESS_H

#include "StochasticProcess.h"
#include "../variables/continuous/NormalRand.h"

/**
 * @brief The WienerProcess class
 */
class RANDLIBSHARED_EXPORT WienerProcess : public StochasticProcess
{
    double mu, var;
public:
    WienerProcess(double mean = 0, double variance = 1);
    void setMean(double mean);
    void setVar(double variance);

    inline double getMean() { return mu; }
    inline double getVar() { return var; }

    bool generate(const QVector<double> &time, QVector<double> &output);

    virtual void E(const QVector<double> &time, QVector<double> &output) const;
    virtual void Var(const QVector<double> &time, QVector<double> &output) const;
};

#endif // WIENERPROCESS_H

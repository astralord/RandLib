#ifndef WIENERPROCESS_H
#define WIENERPROCESS_H

#include "StochasticProcess.h"
#include "../variables/continuous/NormalRand.h"

/**
 * @brief The WienerProcess class
 */
class RANDLIBSHARED_EXPORT WienerProcess : public StochasticProcess
{
public:
    WienerProcess();

    static bool generate(const QVector<double> &time, QVector<double> &output);

    virtual void M(const QVector<double> &time, QVector<double> &output) const;
    virtual void Var(const QVector<double> &time, QVector<double> &output) const;
};

#endif // WIENERPROCESS_H

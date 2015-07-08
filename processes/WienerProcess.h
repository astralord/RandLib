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
    static bool generate(double T, QVector<double> &output);
};

#endif // WIENERPROCESS_H

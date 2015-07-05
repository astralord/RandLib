#ifndef WIENERPROCESS_H
#define WIENERPROCESS_H

#include "StochasticProcess.h"
#include "../variables/continuous/NormalRand.h"

class RANDLIBSHARED_EXPORT WienerProcess : public StochasticProcess
{
public:
    WienerProcess();

    bool generate(const QVector<double> &time, QVector<double> &output);
    bool generate(double T, QVector<double> &output);
};

#endif // WIENERPROCESS_H

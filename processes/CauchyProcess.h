#ifndef CAUCHYPROCESS_H
#define CAUCHYPROCESS_H

#include "StableProcess.h"

class RANDLIBSHARED_EXPORT CauchyProcess : public StableProcess
{
public:
    CauchyProcess(double drift, double volatility, double deltaT = 1.0);
    void nextImpl() override;

    double QuantileImpl(double t, double p) const override;
};


#endif // CAUCHYPROCESS_H

#ifndef BROWNIANMOTION_H
#define BROWNIANMOTION_H

#include "StableProcess.h"
#include "../distributions/univariate/continuous/NormalRand.h"

/**
 * @brief The BrownianMotion class
 */
class RANDLIBSHARED_EXPORT BrownianMotion : public StableProcess
{
    double sqrtDt;
public:
    explicit BrownianMotion(double deltaT = 1.0);

private:
    void nextImpl() override;
    void nextImpl(double deltaT) override;

    double MeanImpl(double t) const override;
    double VarianceImpl(double t) const override;
    double QuantileImpl(double t, double p) const override;
};


#endif // BROWNIANMOTION_H

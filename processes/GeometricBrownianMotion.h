#ifndef GEOMETRICBROWNIANMOTION_H
#define GEOMETRICBROWNIANMOTION_H

#include "StochasticProcess.h"
#include "WienerProcess.h"

class GeometricBrownianMotion : public StochasticProcess
{
    double mu, sigma, S0;
    WienerProcess W;

public:
    GeometricBrownianMotion(double deltaT = 1.0, double drift = 0.0, double volatility = 1.0, double initialValue = 1.0);

    void setParameters(double drift, double volatility, double initialValue);
    inline double getDrift() { return mu; }
    inline double getVolatility() { return sigma; }

    double next() const;
    double next(double deltaT) const;

    void Mean(const QVector<double> &time, QVector<double> &output) const;
    void Variance(const QVector<double> &time, QVector<double> &output) const;
};

#endif // GEOMETRICBROWNIANMOTION_H

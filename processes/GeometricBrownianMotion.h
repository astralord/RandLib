#ifndef GEOMETRICBROWNIANMOTION_H
#define GEOMETRICBROWNIANMOTION_H

#include "StochasticProcess.h"
#include "WienerProcess.h"

class GeometricBrownianMotion : public StochasticProcess
{
    double mu, sigma, S0;
    double generateCoef; /// mu - sigma^2 / 2
public:
    GeometricBrownianMotion(double drift, double volatility, double initialValue);
    void setParameters(double drift, double volatility, double initialValue);
    inline double getDrift() { return mu; }
    inline double getVolatility() { return sigma; }

    bool generate(const QVector<double> &time, QVector<double> &output);

    virtual void M(const QVector<double> &time, QVector<double> &output) const;
    virtual void Var(const QVector<double> &time, QVector<double> &output) const;
};

#endif // GEOMETRICBROWNIANMOTION_H

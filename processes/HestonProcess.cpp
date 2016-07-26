#include "HestonProcess.h"

HestonProcess::HestonProcess(double drift, double volatilityMean, double reversionSpeed, double volatility, double initialValue, double volatilityInitialValue, double correlation, double deltaT) :
    StochasticProcess(deltaT, initialValue),
    mu(drift),
    rho(std::max(0.0, std::min(correlation, 1.0))),
    V(volatilityMean, reversionSpeed, volatility, volatilityInitialValue),
    BW(0, 0, dt, dt, rho)
{

}

void HestonProcess::nextImpl()
{
    /// Euler discretization
    /// Simulate CIR process
    static double v = V.getCurrentValue();
    double sqrtV = std::sqrt(v);
    DoublePair X = BW.variate();
    double dB = X.first, dW  = X.second;
    v += V.getReversionSpeed() * (V.getMean() - v) * dt;
    v += sqrtV * V.getVolatility() * dW;

    /// Simulate Heston from CIR
    double s = 1.0 + mu * dt + sqrtV * dB;
    s = std::max(s, 0.0);
    currentValue *= s;
}

double HestonProcess::MeanImpl(double t) const
{
    // TODO:
    return t;
}

double HestonProcess::VarianceImpl(double t) const
{
    // TODO:
    return t;
}

#include "GeometricBrownianMotion.h"

GeometricBrownianMotion::GeometricBrownianMotion(double drift, double volatility, double initialValue)
{
    setParameters(drift, volatility, initialValue);
}

void GeometricBrownianMotion::setParameters(double drift, double volatility, double initialValue)
{
    mu = drift;
    sigma = std::max(volatility, MIN_POSITIVE);
    S0 = initialValue;
    generateCoef = mu - .5 * sigma * sigma;
}

bool GeometricBrownianMotion::generate(const QVector<double> &time, QVector<double> &output)
{
    int size = std::min(time.size(), output.size());
    if (size <= 0)
        return false;
    if (!WienerProcess::generate(time, output))
        return false;
    // TODO: add all coefs to WienerProcess
    for (int i = 1; i < size; ++i)
    {
        output[i] *= sigma;
        output[i] += generateCoef * time[i];
        output[i] = S0 * std::exp(output[i]);
    }

    return true;
}

void GeometricBrownianMotion::E(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = S0 * std::exp(mu * time[i]);
}

void GeometricBrownianMotion::Var(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
    {
        output[i] = S0 * std::exp(mu * time[i]);
        output[i] *= output[i];
        output[i] *= std::expm1(sigma * sigma * time[i]);
    }
}

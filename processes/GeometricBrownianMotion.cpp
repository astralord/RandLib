#include "GeometricBrownianMotion.h"

GeometricBrownianMotion::GeometricBrownianMotion(double drift, double volatility, double initialValue)
{
    setParameters(drift, volatility, initialValue);
}

void GeometricBrownianMotion::setParameters(double drift, double volatility, double initialValue)
{
    mu = drift;
    sigma = volatility;
    S0 = initialValue;
    generateCoef = mu - .5 * .5 * sigma;
}

bool GeometricBrownianMotion::generate(const QVector<double> &time, QVector<double> &output)
{
    int size = std::min(time.size(), output.size());
    if (size <= 0)
        return false;
    WienerProcess::generate(time, output);
    for (int i = 1; i < size; ++i)
    {
        output[i] *= sigma;
        output[i] += generateCoef * time[i];
        output[i] = S0 * std::exp(output[i]);
    }

    return true;
}

bool GeometricBrownianMotion::generate(double T, QVector<double> &output)
{
    int size = output.size();
    if (size <= 0)
        return false;
    WienerProcess::generate(T, output);
    if (size <= 1)
        return true;
    double deltaT = T / (size - 1);
    for (int i = 1; i < size; ++i)
    {
        output[i] *= sigma;
        output[i] += generateCoef * deltaT * i;
        output[i] = S0 * std::exp(output[i]);
    }

    return true;
}

#include "WienerProcess.h"

WienerProcess::WienerProcess(double deltaT, double mean, double variance) :
    StochasticProcess(deltaT)
{
    setMean(mean);
    setVar(variance);
}

void WienerProcess::setMean(double mean)
{
    mu = mean;
    muT = mu * dt;
}

void WienerProcess::setVar(double variance)
{
    var = variance;
    if (var <= 0)
        var = 1.0;
    sigma = std::sqrt(var);
    sigmaT = std::sqrt(var * dt);
}

void WienerProcess::setSigma(double volatility)
{
    sigma = volatility;
    if (sigma <= 0)
        sigma = 1.0;
    var = sigma * sigma;
    sigmaT = sigma * std::sqrt(dt);
}

double WienerProcess::next() const
{
    lastValue += NormalRand::variate(muT, sigmaT);
    return lastValue;
}

double WienerProcess::next(double deltaT) const
{
    lastValue += NormalRand::variate(mu * deltaT, std::sqrt(var * deltaT));
    return lastValue;
}

void WienerProcess::Mean(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = mu * time[i];
}

void WienerProcess::Variance(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = time[i] * var;
}

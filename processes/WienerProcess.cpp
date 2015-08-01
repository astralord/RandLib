#include "WienerProcess.h"

WienerProcess::WienerProcess(double mean, double variance)
{
    setMean(mean);
    setVar(variance);
}

void WienerProcess::setMean(double mean)
{
    mu = mean;
}

void WienerProcess::setVar(double variance)
{
    var = std::sqrt(std::max(variance, MIN_POSITIVE));
}

bool WienerProcess::generate(const QVector<double> &time, QVector<double> &output)
{
    int size = std::min(time.size(), output.size());
    if (size <= 0)
        return false;
    output[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        double deltaT = time[i] - time[i - 1];
        output[i] = NormalRand::variate(output[i - 1] + mu * deltaT, std::sqrt(var * deltaT));
    }
    return true;
}

void WienerProcess::E(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = mu * time[i];
}

void WienerProcess::Var(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = time[i] * var;
}

#include "WienerProcess.h"

WienerProcess::WienerProcess()
{
}

bool WienerProcess::generate(const QVector<double> &time, QVector<double> &output)
{
    int size = std::min(time.size(), output.size());
    if (size <= 0)
        return false;
    output[0] = 0;
    NormalRand rv;
    for (int i = 1; i < size; ++i)
    {
        rv.setVar(time[i] - time[i - 1]);
        output[i] = output[i - 1] + rv.value();
    }
    return true;
}

void WienerProcess::M(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = 0;
}

void WienerProcess::Var(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = time[i];
}

#include "RandomVariable.h"

RandomVariable::RandomVariable()
{
}

void RandomVariable::sample(QVector<double> &outputData)
{
    for (double &var : outputData)
        var = variate();
}

void RandomVariable::cdf(const QVector<double> &x, QVector<double> &y)
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

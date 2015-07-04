#include "RandomVariable.h"

template< typename T>
RandomVariable<T>::RandomVariable()
{
}

template< typename T>
void RandomVariable<T>::sample(QVector<T> &outputData)
{
    for (T &var : outputData)
        var = value();
}

template< typename T>
void RandomVariable<T>::cdf(const QVector<double> &x, QVector<double> &y)
{
    if (x.size() != y.size())
        return;
    for (int i = 0; i != x.size(); ++i)
        y[i] = F(x[i]);
}

template class RandomVariable<int>;
template class RandomVariable<double>;

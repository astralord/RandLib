#include "RandomVariable.h"

template< typename T>
RandomVariable<T>::RandomVariable()
{
}

template< typename T>
void RandomVariable<T>::sample(QVector<T> &outputData)
{
    for (T &var : outputData)
        var = variate();
}

template< typename T>
void RandomVariable<T>::cdf(const QVector<double> &x, QVector<double> &y)
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

template class RandomVariable<int>;
template class RandomVariable<double>;

#include "DiscreteRand.h"

template< typename T>
void DiscreteRand<T>::pmf(const QVector<T> &x, QVector<double> &y) const
{
    if (x.size() != y.size())
        return;
    for (int i = 0; i != x.size(); ++i)
        y[i] = P(x[i]);
}

template< typename T>
double DiscreteRand<T>::likelihood(const QVector<T> &sample) const
{
    double res = 0.0;
    for (int i = 0; i != sample.size(); ++i)
        res += std::log(P(sample[i]));
    return res;
}

template class DiscreteRand<int>;
template class DiscreteRand<double>;

#include <RandomVariable.h>
#include <time.h>
#include <QDebug>

template< typename T>
unsigned long RandomVariable<T>::startPoint = 123456789;

template< typename T>
RandomVariable<T>::RandomVariable()
{
    startPoint ^= time(0);
    X = Y = Z = W = startPoint;
    C = 0;
    startPoint = fastKISS();

    N = 0;
    CODE = 1;
    Q = 0;
}

template< typename T>
unsigned long RandomVariable<T>::fastKISS()
{
    Y ^= Y << 5;
    Y ^= Y >> 7;
    Y ^= Y << 22;

    int t = Z + W + C;
    Z = W;
    C = t < 0;
    W = t & 2147483647;
    X += 1411392427;

    return X + Y + W;
}

template< typename T>
unsigned long RandomVariable<T>::quasiGen()
{
    if (N > 0)
    {
        unsigned value = N;
        CODE = 1;
        while (value & 1) {
            value >>= 1;
            ++CODE;
        }

        Q ^= 1 << (32 - CODE);
    }

    ++N;
    return Q;
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

/// CONTINUOUS

void ContinuousRand::pdf(const QVector<double> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

double ContinuousRand::likelihood(const QVector<double> &sample) const
{
    double res = 1.0;
    for (int i = 0; i != sample.size(); ++i)
        res *= f(sample[i]);
    return res;
}

double ContinuousRand::loglikelihood(const QVector<double> &sample) const
{
    double res = 0.0;
    for (int i = 0; i != sample.size(); ++i)
        res += std::log(f(sample[i]));
    return res;
}

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

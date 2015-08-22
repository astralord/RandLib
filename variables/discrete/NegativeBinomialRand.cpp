#include "NegativeBinomialRand.h"

template< typename T>
NegativeBinomialRand<T>::NegativeBinomialRand(T number, double probability) :
    G(probability)
{
    setNumber(number);
    setProbability(probability);
    setName();
}

template< typename T>
void NegativeBinomialRand<T>::setName()
{
    nameStr = "Negative Binomial(" + toStringWithPrecision(getNumber()) + ", " + toStringWithPrecision(getProbability()) + ")";
}

template< typename T>
void NegativeBinomialRand<T>::setProbability(double probability)
{
    r = std::min(std::max(probability, MIN_POSITIVE), 1.0);
    G.setProbability(r);
}

template< typename T>
void NegativeBinomialRand<T>::setNumber(T number)
{
    k = std::max(number, static_cast<T>(1.0));
}

template< typename T>
double NegativeBinomialRand<T>::P(int k) const
{
    return k;
}

template< typename T>
double NegativeBinomialRand<T>::F(double x) const
{
    return x;
}

template< typename T>
double NegativeBinomialRand<T>::variate() const
{
    double res = 0;
    for (int i = 0; i != static_cast<int>(k); ++i)
        res += G.variate();
    return res;
}

template class NegativeBinomialRand<int>;
template class NegativeBinomialRand<double>;

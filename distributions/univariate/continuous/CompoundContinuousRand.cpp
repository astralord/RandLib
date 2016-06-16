#include "CompoundContinuousRand.h"

template <typename T>
double CompoundContinuousRand<T>::f(double x) const
{
    return X.ExpectedValue([this, x](double t)
    {
        return Y.f(x - t);
    }, Mean());
}

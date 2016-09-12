#include "StochasticProcess.h"

template <typename T>
StochasticProcess<T>::StochasticProcess(double deltaT, T initialValue) :
    currentTime(0.0),
    dt(std::max(deltaT, 0.0)),
    currentValue(initialValue)
{
}

template <typename T>
void StochasticProcess<T>::reset(T initialValue)
{
    currentTime = 0.0;
    currentValue = initialValue;
}

template <typename T>
T StochasticProcess<T>::Next()
{
    currentTime += dt;
    NextImpl();
    return currentValue;
}

template <typename T>
double StochasticProcess<T>::Mean(double t) const
{
    if (t < currentTime)
        return NAN; /// we consider only future time
    return (t > currentTime) ? MeanImpl(t) : currentValue;
}

template <typename T>
double StochasticProcess<T>::Variance(double t) const
{
    if (t < currentTime)
        return NAN; /// we consider only future time
    return (t > currentTime) ? VarianceImpl(t) : 0.0;
}

template class StochasticProcess<double>;
template class StochasticProcess<int>;
//template class StochasticProcess<DoublePair>;

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
T StochasticProcess<T>::next()
{
    currentTime += dt;
    nextImpl();
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

template <typename T>
double StochasticProcess<T>::Quantile(double t, double p) const
{
    if (p < 0 || p > 1 || t < currentTime)
        return NAN;
    return QuantileImpl(t, p);
}

template <typename T>
void StochasticProcess<T>::Quantile(const std::vector<double> &t, std::vector<double> &outputData, double p) const
{
    /// Assumed that t is sorted
    if (p < 0 || p > 1 || t[0] < currentTime)
        return;
    if (t.size() > outputData.size())
        return;
    return this->QuantileImpl(t, outputData, p);
}

template <typename T>
void StochasticProcess<T>::QuantileImpl(const std::vector<double> &t, std::vector<double> &outputData, double p) const
{
    for (size_t i = 0; i != t.size(); ++i)
        outputData[i] = QuantileImpl(t[i], p);
}


template class StochasticProcess<double>;
template class StochasticProcess<int>;
//template class StochasticProcess<DoublePair>;

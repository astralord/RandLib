#include "JumpDiffusionProcess.h"

template <typename T>
JumpDiffusionProcess<T>::JumpDiffusionProcess(double drift, double volatility, double rate, const UnivariateProbabilityDistribution<T> &jumpDistribution, double deltaT) :
    StochasticProcess(deltaT),
    B(drift, volatility, dt),
    J(rate, jumpDistribution, dt)
{

}

template <typename T>
void JumpDiffusionProcess<T>::NextImpl()
{
    currentValue = B.Next() + J.Next();
}

template <typename T>
double JumpDiffusionProcess<T>::MeanImpl(double t) const
{
    return B.Mean(t) + J.Mean(t);
}

template <typename T>
double JumpDiffusionProcess<T>::VarianceImpl(double t) const
{
    return B.Variance(t) + J.Variance(t);
}

template class JumpDiffusionProcess<double>;
template class JumpDiffusionProcess<int>;

#include "JumpDiffusionProcess.h"

template <typename T>
JumpDiffusionProcess<T>::JumpDiffusionProcess(double drift, double volatility, double rate, const UnivariateProbabilityDistribution<T> &jumpDistribution, double deltaT) :
    StochasticProcess(deltaT),
    B(drift, volatility, dt),
    J(rate, jumpDistribution, dt)
{

}

template <typename T>
void JumpDiffusionProcess<T>::nextImpl()
{
    currentValue = B.next() + J.next();
}

template <typename T>
void JumpDiffusionProcess<T>::nextImpl(double deltaT)
{
    currentValue = B.next(deltaT) + J.next(deltaT);
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

template <typename T>
double JumpDiffusionProcess<T>::QuantileImpl(double t, double p) const
{
    // very rough approximation by normal distribution
    // can work only for big t and existing mean and variance
    NormalRand X(Mean(t), Variance(t));
    return X.Quantile(p);
}

template class JumpDiffusionProcess<double>;
template class JumpDiffusionProcess<int>;

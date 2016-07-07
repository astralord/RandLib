#include "CompoundPoissonProcess.h"

template <typename T>
CompoundPoissonProcess<T>::CompoundPoissonProcess(double rate, const UnivariateProbabilityDistribution<T> &jumpDistribution, double deltaT) :
    StochasticProcess(deltaT),
    N(rate, deltaT),
    Y(jumpDistribution),
    jumpsAmount(0)
{
}

template <typename T>
void CompoundPoissonProcess<T>::nextImpl()
{
    N.next();
    double Nt = N.getCurrentValue();
    while (Nt > jumpsAmount) {
        currentValue += Y.variate();
        ++jumpsAmount;
    }
}

template <typename T>
void CompoundPoissonProcess<T>::nextImpl(double)
{
    nextImpl();
}

template <typename T>
double CompoundPoissonProcess<T>::MeanImpl(double t) const
{
    return currentValue + Y.Mean() * N.getRate() * (t - currentTime);
}

template <typename T>
double CompoundPoissonProcess<T>::VarianceImpl(double t) const
{
    return Y.Variance() * N.getRate() * (t - currentTime);
}

template <typename T>
double CompoundPoissonProcess<T>::QuantileImpl(double t, double p) const
{
    // TODO
    // Quantile is very hard to derive - I'll leave it for now
    return Mean(t) + (p - 0.5) * 6 * Variance(t);
}


template class CompoundPoissonProcess<double>;
template class CompoundPoissonProcess<int>;

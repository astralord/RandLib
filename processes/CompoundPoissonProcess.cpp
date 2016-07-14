#include "CompoundPoissonProcess.h"
#include "../distributions/univariate/continuous/NormalRand.h"

template <typename T>
CompoundPoissonProcess<T>::CompoundPoissonProcess(double rate, const UnivariateProbabilityDistribution<T> &jumpDistribution, double deltaT) :
    StochasticProcess<T>(deltaT),
    N(rate, this->dt),
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
        this->currentValue += Y.variate();
        ++jumpsAmount;
    }
}

template <typename T>
double CompoundPoissonProcess<T>::MeanImpl(double t) const
{
    return this->currentValue + Y.Mean() * N.getRate() * (t - this->currentTime);
}

template <typename T>
double CompoundPoissonProcess<T>::VarianceImpl(double t) const
{
    return Y.Variance() * N.getRate() * (t - this->currentTime);
}

template class CompoundPoissonProcess<double>;
template class CompoundPoissonProcess<int>;

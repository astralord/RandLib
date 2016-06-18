#include "SumRand.h"

template < typename T >
SumContinuousRand<T>::SumContinuousRand(const ContinuousDistribution &leftRV, const UnivariateProbabilityDistribution<T> &rightRV)
    : X(leftRV), Y(rightRV)
{
}

template < typename T >
std::string SumContinuousRand<T>::name() const
{
    return X.name() + Y.name();
}

template < typename T >
double SumContinuousRand<T>::f(double x) const
{
    double startPoint = Mean();
    if (!std::isfinite(startPoint))
        startPoint = 0.0;
    return Y.ExpectedValue([this, x](double t) {
        return X.f(x - t);
    }, startPoint);
}

template < typename T >
double SumContinuousRand<T>::F(double x) const
{
    double startPoint = Mean();
    if (!std::isfinite(startPoint))
        startPoint = 0.0;
    return Y.ExpectedValue([this, x](double t)
    {
        return X.F(x - t);
    }, startPoint);
}

template < typename T >
double SumContinuousRand<T>::variate() const
{
    return X.variate() + Y.variate();
}

template < typename T >
double SumContinuousRand<T>::Mean() const
{
    return X.Mean() + Y.Mean();
}

template < typename T >
double SumContinuousRand<T>::Variance() const
{
    return X.Variance() + Y.Variance();
}

template < typename T >
std::complex<double> SumContinuousRand<T>::CF(double t) const
{
    return X.CF(t) * Y.CF(t);
}

template < typename T >
double SumContinuousRand<T>::Skewness() const
{
    double skewX = X.Skewness();
    double skewY = Y.Skewness();
    double varX = X.Variance();
    double varY = Y.Variance();
    double skew = skewX * std::pow(varX, 1.5);
    skew += skewY * std::pow(varY, 1.5);
    return skew / std::pow(varX + varY, 1.5);
}

template < typename T >
double SumContinuousRand<T>::ExcessKurtosis() const
{
    double varXSq = X.Variance();
    varXSq *= varXSq;
    double varYSq = Y.Variance();
    varYSq *= varYSq;
    double kurtX = X.Kurtosis();
    double kurtY = Y.Kurtosis();
    double kurt = kurtX * varXSq * varXSq;
    kurt += kurtY * varYSq * varYSq;
    kurt += 6 * varXSq * varYSq;
    double denom = varXSq + varYSq;
    return kurt / (denom * denom);
}

template class SumContinuousRand<double>;
template class SumContinuousRand<int>;

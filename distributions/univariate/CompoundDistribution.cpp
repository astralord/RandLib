#include "CompoundDistribution.h"

template < typename T1, typename T2 >
CompoundDistribution<T1, T2>::CompoundDistribution(const UnivariateProbabilityDistribution<T1> &leftRV, const UnivariateProbabilityDistribution<T2> &rightRV)
    : X(leftRV), Y(rightRV)
{
}

template < typename T1, typename T2 >
T1 CompoundDistribution<T1, T2>::variate() const
{
    return X.variate() + Y.variate();
}

template < typename T1, typename T2 >
double CompoundDistribution<T1, T2>::Mean() const
{
    return X.Mean() + Y.Mean();
}

template < typename T1, typename T2 >
double CompoundDistribution<T1, T2>::Variance() const
{
    return X.Variance() + Y.Variance();
}

template < typename T1, typename T2 >
std::complex<double> CompoundDistribution<T1, T2>::CF(double t) const
{
    return X.CF(t) * Y.CF(t);
}

template < typename T1, typename T2 >
double CompoundDistribution<T1, T2>::Skewness() const
{
    double skewX = X.Skewness();
    double skewY = Y.Skewness();
    double varX = X.Variance();
    double varY = Y.Variance();
    double skew = skewX * std::pow(varX, 1.5);
    skew += skewY * std::pow(varY, 1.5);
    return skew / std::pow(varX + varY, 1.5);
}

template < typename T1, typename T2 >
double CompoundDistribution<T1, T2>::ExcessKurtosis() const
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

template class CompoundRand<double, double>;
template class CompoundRand<double, int>;
template class CompoundRand<int, int>;

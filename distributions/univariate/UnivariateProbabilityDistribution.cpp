#include "UnivariateProbabilityDistribution.h"

template< typename T >
UnivariateProbabilityDistribution<T>::UnivariateProbabilityDistribution()
{
}

template< typename T >
void UnivariateProbabilityDistribution<T>::QuantileFunction(const std::vector<double> &p, std::vector<double> &y)
{
    int size = std::min(p.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = Quantile(p[i]);
}

template< typename T >
std::complex<double> UnivariateProbabilityDistribution<T>::CF(double t) const
{
    if (t == 0)
        return std::complex<double>(1, 0);
    double startPoint = Mean();
    if (!std::isfinite(startPoint))
        startPoint = Median();
    double re = ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x);
    }, startPoint);
    double im = ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    }, startPoint);

    return std::complex<double>(re, im);
}

template< typename T >
void UnivariateProbabilityDistribution<T>::CharacteristicFunction(const std::vector<double> &t, std::vector<std::complex<double> > &y) const
{
    int size = std::min(t.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = CF(t[i]);
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Median() const
{
    return Quantile(0.5);
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Skewness() const
{
    double mu = Mean();
    if (!std::isfinite(mu))
        return NAN;

    double var = Variance();
    if (!std::isfinite(var))
        return NAN;

    double sum = ExpectedValue([this, mu] (double x)
    {
        double skew = x - mu;
        return skew * skew * skew;
    }, mu);

    return sum / (var * std::sqrt(var));
}

template< typename T >
double UnivariateProbabilityDistribution<T>::ExcessKurtosis() const
{
    double mu = Mean();
    if (!std::isfinite(mu))
        return NAN;

    double var = Variance();
    if (!std::isfinite(var))
        return NAN;

    double sum = ExpectedValue([this, mu] (double x)
    {
        double kurtosis = x - mu;
        kurtosis *= kurtosis;
        return kurtosis * kurtosis;
    }, mu);

    return sum / (var * var) - 3;
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Kurtosis() const
{
    return ExcessKurtosis() + 3.0;
}

template class UnivariateProbabilityDistribution<double>;
template class UnivariateProbabilityDistribution<int>;

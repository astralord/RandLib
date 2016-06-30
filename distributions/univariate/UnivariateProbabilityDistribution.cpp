#include "UnivariateProbabilityDistribution.h"

template< typename T >
UnivariateProbabilityDistribution<T>::UnivariateProbabilityDistribution()
{
}

template< typename T >
bool UnivariateProbabilityDistribution<T>::isLeftBounded() const
{
    SUPPORT_TYPE supp = supportType();
    return (supp == RIGHTSEMIFINITE_T || supp == FINITE_T);
}

template< typename T >
bool UnivariateProbabilityDistribution<T>::isRightBounded() const
{
    SUPPORT_TYPE supp = supportType();
    return (supp == LEFTSEMIFINITE_T || supp == FINITE_T);
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

    if (std::fabs(t) < 0.0001)
    {
        double mean4 = FourthMoment();
        if (std::isfinite(mean4)) {
            double mean3 = ThirdMoment();
            double mean2 = SecondMoment();
            double mean1 = Mean();

            double tSq = t * t;
            double re = 1 - 0.5 * mean2 * tSq + mean4 * tSq * tSq / 24;
            double im = mean1 - mean3 * tSq / 6;
            im *= t;
            return std::complex<double>(re, im);
        }
    }

    double startPoint = Mean();
    if (!std::isfinite(startPoint)) {
        if (isLeftBounded())
            startPoint = MinValue();
        else if (isRightBounded())
            startPoint = MaxValue();
        else {
            startPoint = Median();
            if (!std::isfinite(startPoint))
                startPoint = 0.0;
        }
    }
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
void UnivariateProbabilityDistribution<T>::HazardFunction(const std::vector<double> &x, std::vector<double> &y) const
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = Hazard(x[i]);
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
double UnivariateProbabilityDistribution<T>::SecondMoment() const
{
    double mean = Mean();
    return mean * mean + Variance();
}

template< typename T >
double UnivariateProbabilityDistribution<T>::ThirdMoment() const
{
    double mean = Mean();
    double variance = Variance();
    double skewness = Skewness();

    double moment = skewness * std::sqrt(variance) * variance;
    moment += mean * mean * mean;
    moment += mean * variance;
    return moment;
}

template< typename T >
double UnivariateProbabilityDistribution<T>::FourthMoment() const
{
    double mean = Mean();
    double variance = Variance();
    double moment3 = ThirdMoment();
    double kurtosis = Kurtosis();
    double meanSq = mean * mean;

    double moment = kurtosis * variance * variance;
    moment -= 6 * meanSq * variance;
    moment -= 3 * meanSq * meanSq;
    moment += 4 * mean * moment3;
    return moment;
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Kurtosis() const
{
    return ExcessKurtosis() + 3.0;
}

template class UnivariateProbabilityDistribution<double>;
template class UnivariateProbabilityDistribution<int>;

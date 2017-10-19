#include "UnivariateProbabilityDistribution.h"

template< typename T >
UnivariateProbabilityDistribution<T>::UnivariateProbabilityDistribution()
{
}

template< typename T >
bool UnivariateProbabilityDistribution<T>::isLeftBounded() const
{
    SUPPORT_TYPE supp = SupportType();
    return (supp == RIGHTSEMIFINITE_T || supp == FINITE_T);
}

template< typename T >
bool UnivariateProbabilityDistribution<T>::isRightBounded() const
{
    SUPPORT_TYPE supp = SupportType();
    return (supp == LEFTSEMIFINITE_T || supp == FINITE_T);
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Quantile(double p) const
{
    if (p < 0.0 || p > 1.0)
        return NAN;
    double minVal = this->MinValue();
    if (p == 0.0)
        return minVal;
    double maxVal = this->MaxValue();
    if (p == 1.0)
        return maxVal;
    double x = this->quantileImpl(p);
    if (x < minVal)
        return minVal;
    return (x > maxVal) ? maxVal : x;
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Quantile1m(double p) const
{
    if (p < 0.0 || p > 1.0)
        return NAN;
    double minVal = this->MinValue();
    if (p == 1.0)
        return minVal;
    double maxVal = this->MaxValue();
    if (p == 0.0)
        return maxVal;
    double x = this->quantileImpl1m(p);
    if (x < minVal)
        return minVal;
    return (x > maxVal) ? maxVal : x;
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
    if (t == 0.0)
        return 1.0;
    return (t > 0.0) ? this->CFImpl(t) : std::conj(this->CFImpl(-t));
}

template< typename T >
std::complex<double> UnivariateProbabilityDistribution<T>::CFImpl(double t) const
{
    T leftBound = this->MinValue(), rightBound = this->MaxValue();
    if (leftBound == rightBound)
        return std::complex<double>(std::cos(t * leftBound), std::sin(t * leftBound));

    double re = ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x);
    }, leftBound, rightBound);

    double im = ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    }, leftBound, rightBound);

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
T UnivariateProbabilityDistribution<T>::Median() const
{
    return quantileImpl(0.5);
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Skewness() const
{
    double var = Variance();
    if (!std::isfinite(var))
        return NAN;
    double mu = Mean(); /// var is finite, so is mu

    double sum = ExpectedValue([this, mu] (double x)
    {
        double xmmu = x - mu;
        double skewness = xmmu * xmmu * xmmu;
        return skewness;
    }, this->MinValue(), this->MaxValue());

    return sum / std::pow(var, 1.5);
}

template< typename T >
double UnivariateProbabilityDistribution<T>::ExcessKurtosis() const
{
    double var = Variance();
    if (!std::isfinite(var))
        return NAN;
    double mu = Mean(); /// var is finite, so is mu

    double sum = ExpectedValue([this, mu] (double x)
    {
        double xmmu = x - mu;
        double kurtosisSqrt = xmmu * xmmu;
        double kurtosis = kurtosisSqrt * kurtosisSqrt;
        return kurtosis;
    }, this->MinValue(), this->MaxValue());

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
    moment += 3 * mean * variance;
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
bool UnivariateProbabilityDistribution<T>::allElementsAreNotBiggerThan(T value, const std::vector<T> &sample)
{
    for (const T & var : sample) {
        if (var > value)
            return false;
    }
    return true;
}

template< typename T >
bool UnivariateProbabilityDistribution<T>::allElementsAreNotLessThan(T value, const std::vector<T> &sample)
{
    for (const T & var : sample) {
        if (var < value)
            return false;
    }
    return true;
}

template< typename T >
bool UnivariateProbabilityDistribution<T>::allElementsAreNonNegative(const std::vector<T> &sample)
{
    return allElementsAreNotLessThan(0, sample);
}

template< typename T >
bool UnivariateProbabilityDistribution<T>::allElementsArePositive(const std::vector<T> &sample)
{
    for (const T & var : sample) {
        if (var <= 0)
            return false;
    }
    return true;
}

template< typename T >
double UnivariateProbabilityDistribution<T>::Kurtosis() const
{
    return ExcessKurtosis() + 3.0;
}

template< typename T >
double UnivariateProbabilityDistribution<T>::GetSampleSum(const std::vector<T> &sample)
{
    return std::accumulate(sample.begin(), sample.end(), 0.0);
}

template< typename T >
double UnivariateProbabilityDistribution<T>::GetSampleMean(const std::vector<T> &sample)
{
    size_t n = sample.size();
    return (n > 0) ? GetSampleSum(sample) / n : 0.0;
}

template< typename T >
double UnivariateProbabilityDistribution<T>::GetSampleVariance(const std::vector<T> &sample, double mean)
{
    long double sum = 0.0l;
    for (const T & var : sample) {
        double temp = var - mean;
        sum += temp * temp;
    }
    return sum / sample.size();
}

template< typename T >
DoublePair UnivariateProbabilityDistribution<T>::GetSampleMeanAndVariance(const std::vector<T> &sample)
{
    /// Welford's stable method
    long double m = 0.0l, v = 0.0l;
    size_t n = sample.size();
    for (size_t i = 0; i < n; ++i) {
        long double diff = sample[i] - m;
        m += diff / (i + 1);
        v += diff * (sample[i] - m);
    }
    return std::make_pair(m, v / n);
}

template< typename T >
std::tuple<double, double, double, double> UnivariateProbabilityDistribution<T>::GetSampleStatistics(const std::vector<T> &sample)
{
    // TODO: implement method of Knuth and Welford for skewness and kurtosis also!
    DoublePair mv = GetSampleMeanAndVariance(sample);
    size_t n = sample.size();
    long double skewness = 0.0l, kurtosis = 0.0l;
    for (double var : sample) {
        skewness += std::pow(var - mv.first, 3);
        kurtosis += std::pow(var - mv.first, 4);
    }
    skewness /= (n * std::pow(mv.second, 1.5));
    kurtosis /= (n * std::pow(mv.second, 2));
    return std::make_tuple(mv.first, mv.second, skewness, kurtosis);
}

template class UnivariateProbabilityDistribution<double>;
template class UnivariateProbabilityDistribution<int>;

#include "UnivariateDistribution.h"

template< typename T >
UnivariateDistribution<T>::UnivariateDistribution()
{
}

template< typename T >
bool UnivariateDistribution<T>::isLeftBounded() const
{
    SUPPORT_TYPE supp = SupportType();
    return (supp == RIGHTSEMIFINITE_T || supp == FINITE_T);
}

template< typename T >
bool UnivariateDistribution<T>::isRightBounded() const
{
    SUPPORT_TYPE supp = SupportType();
    return (supp == LEFTSEMIFINITE_T || supp == FINITE_T);
}

template< typename T >
double UnivariateDistribution<T>::Quantile(double p) const
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
double UnivariateDistribution<T>::Quantile1m(double p) const
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
void UnivariateDistribution<T>::QuantileFunction(const std::vector<double> &p, std::vector<double> &y)
{
    int size = std::min(p.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = Quantile(p[i]);
}

template< typename T >
std::complex<double> UnivariateDistribution<T>::CF(double t) const
{
    if (t == 0.0)
        return 1.0;
    return (t > 0.0) ? this->CFImpl(t) : std::conj(this->CFImpl(-t));
}

template< typename T >
std::complex<double> UnivariateDistribution<T>::CFImpl(double t) const
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
void UnivariateDistribution<T>::CharacteristicFunction(const std::vector<double> &t, std::vector<std::complex<double> > &y) const
{
    int size = std::min(t.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = CF(t[i]);
}

template< typename T >
void UnivariateDistribution<T>::HazardFunction(const std::vector<double> &x, std::vector<double> &y) const
{
    int size = std::min(x.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = Hazard(x[i]);
}

template< typename T >
T UnivariateDistribution<T>::Median() const
{
    return quantileImpl(0.5);
}

template< typename T >
double UnivariateDistribution<T>::Skewness() const
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
double UnivariateDistribution<T>::ExcessKurtosis() const
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
double UnivariateDistribution<T>::SecondMoment() const
{
    double mean = Mean();
    return mean * mean + Variance();
}

template< typename T >
double UnivariateDistribution<T>::ThirdMoment() const
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
double UnivariateDistribution<T>::FourthMoment() const
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
bool UnivariateDistribution<T>::allElementsAreNotBiggerThan(T value, const std::vector<T> &sample)
{
    for (const T & var : sample) {
        if (var > value)
            return false;
    }
    return true;
}

template< typename T >
bool UnivariateDistribution<T>::allElementsAreNotLessThan(T value, const std::vector<T> &sample)
{
    for (const T & var : sample) {
        if (var < value)
            return false;
    }
    return true;
}

template< typename T >
bool UnivariateDistribution<T>::allElementsAreNonNegative(const std::vector<T> &sample)
{
    return allElementsAreNotLessThan(0, sample);
}

template< typename T >
bool UnivariateDistribution<T>::allElementsArePositive(const std::vector<T> &sample)
{
    for (const T & var : sample) {
        if (var <= 0)
            return false;
    }
    return true;
}

template< typename T >
double UnivariateDistribution<T>::Kurtosis() const
{
    return ExcessKurtosis() + 3.0;
}

template< typename T >
double UnivariateDistribution<T>::GetSampleSum(const std::vector<T> &sample)
{
    return std::accumulate(sample.begin(), sample.end(), 0.0);
}

template< typename T >
double UnivariateDistribution<T>::GetSampleMean(const std::vector<T> &sample)
{
    size_t n = sample.size();
    return (n > 0) ? GetSampleSum(sample) / n : 0.0;
}

template< typename T >
double UnivariateDistribution<T>::GetSampleVariance(const std::vector<T> &sample, double mean)
{
    long double sum = 0.0l;
    for (const T & var : sample) {
        double temp = var - mean;
        sum += temp * temp;
    }
    return sum / sample.size();
}

template< typename T >
DoublePair UnivariateDistribution<T>::GetSampleMeanAndVariance(const std::vector<T> &sample)
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
std::tuple<double, double, double, double> UnivariateDistribution<T>::GetSampleStatistics(const std::vector<T> &sample)
{
    /// Terriberry's extension for skewness and kurtosis
    long double M1{}, M2{}, M3{}, M4{};
    long double m1{}, m2{}, m3{}, m4{};
    size_t n = sample.size(), k = 0;
    size_t t = 0;
    static constexpr size_t BIG_NUMBER = 10000;
    for (const T & var : sample) {
        ++k;
        long double delta = var - m1;
        long double delta_k = delta / k;
        long double delta_kSq = delta_k * delta_k;
        long double term1 = delta * delta_k * (k - 1);
        m1 += delta_k;
        m4 += term1 * delta_kSq * (k * k - 3 * k + 3) + 6 * delta_kSq * m2 - 4 * delta_k * m3;
        m3 += term1 * delta_k * (k - 2) - 3 * delta_k * m2;
        m2 += term1;

        /// This looks like a hack and unfortunately it is. The reason of it is that the algorithm
        /// can become unstable for sufficiently large k. For now we restart the algorithm when
        /// k reaches some big number. In the future this can be parallelized for faster implementation.
        if (k >= BIG_NUMBER) {
            long double Delta = m1 - M1;
            long double DeltaSq = Delta * Delta;
            size_t tp1 = t + 1;
            /// M4
            M4 += m4 + (DeltaSq * DeltaSq * t * (t * t - t + 1) * BIG_NUMBER) / (tp1 * tp1 * tp1);
            M4 += 6 * DeltaSq * (t * t * m2 + M2) / (tp1 * tp1);
            M4 += 4 * Delta * (t * m3 - M3) / tp1;
            /// M3
            M3 += m3 + (DeltaSq * Delta * t * (t - 1) * BIG_NUMBER) / (tp1 * tp1);
            M3 += 3 * Delta * (t * m2 - M2) / tp1;
            /// M2 and M1
            M2 += m2 + (DeltaSq * t) / tp1;
            M1 += Delta / tp1;
            k = 0;
            m1 = 0;
            m2 = 0;
            m3 = 0;
            m4 = 0;
            ++t;
        }
    }

    /// If something left - add the residue
    double res = static_cast<double>(k) / BIG_NUMBER;
    if (res != 0) {
        long double Delta = m1 - M1;
        long double DeltaSq = Delta * Delta;
        double tpres = t + res;
        /// M4
        M4 += m4 + (DeltaSq * DeltaSq * t * res * (t * t - res * t + res * res) * BIG_NUMBER) / (tpres * tpres * tpres);
        M4 += 6 * DeltaSq * (t * t * m2 + res * res * M2) / (tpres * tpres);
        M4 += 4 * Delta * (t * m3 - res * M3) / tpres;
        /// M3
        M3 += m3 + (DeltaSq * Delta * t * res * (t - res) * BIG_NUMBER) / (tpres * tpres);
        M3 += 3 * Delta * (t * m2 - res * M2) / tpres;
        /// M2 and M1
        M2 += m2 + (DeltaSq * t * res) / tpres;
        M1 += Delta * res / tpres;
    }

    double variance = M2 / n;
    double skewness = std::sqrt(n) * M3 / std::pow(M2, 1.5);
    double exkurtosis = (n * M4) / (M2 * M2) - 3.0;
    return std::make_tuple(M1, variance, skewness, exkurtosis);
}

template class UnivariateDistribution<double>;
template class UnivariateDistribution<int>;

#pragma once

#include <span>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/distribution.hpp>

namespace randlib
{

/// Shape of the support of a univariate distribution
enum class Support
{
    Finite, ///< [a, b]
    RightSemiFinite, ///< [a, ∞)
    LeftSemiFinite, ///< (-∞, b]
    Infinite ///< (-∞, ∞)
};

/**
 * @brief The UnivariateDistribution class <BR>
 * Abstract class for all univariate probability distributions
 */
template <typename T>
class UnivariateDistribution : public ProbabilityDistribution<T>
{
protected:
    UnivariateDistribution();

public:
    virtual ~UnivariateDistribution() {}

protected:
public:
    /**
     * @return type of support
     */
    [[nodiscard]] virtual Support SupportType() const = 0;

    /**
     * @return true if distribution is bounded from the left
     */
    [[nodiscard]] bool isLeftBounded() const;

    /**
     * @return true if distribution is bounded from the right
     */
    [[nodiscard]] bool isRightBounded() const;

    /**
     * @return Mathematical expectation
     */
    [[nodiscard]] virtual long double Mean() const = 0;

    /**
     * @return Variance of random variable
     */
    [[nodiscard]] virtual long double Variance() const = 0;

private:
    /**
     * @param p
     * @param initValue initial value of x
     * @return such x that F(x) = p
     */
    virtual T quantileImpl(double p, T initValue) const = 0;
    virtual T quantileImpl(double p) const = 0;

    /**
     * @param p
     * @param initValue initial value of x
     * @return such x that F(x) = 1 - p
     */
    virtual T quantileImpl1m(double p, T initValue) const = 0;
    virtual T quantileImpl1m(double p) const = 0;

protected:
    /**
     * @param t
     * @return characteristic function implementation (inverse Fourier transform of pdf)
     */
    virtual std::complex<double> CFImpl(double t) const;

    /**
     * @param funPtr pointer on function g(x) with finite support which expected value should be returned
     * @param minPoint min{x | g(x) ≠ 0}
     * @param maxPoint max{x | g(x) ≠ 0}
     * @return E[g(x)]
     */
    virtual long double ExpectedValue(const std::function<double(T)> &funPtr, T minPoint, T maxPoint) const = 0;

public:
    /**
     * @param p
     * @return quantileImpl(p) if p is in (0, 1)
     */
    [[nodiscard]] T Quantile(double p) const;

    /**
     * @param p
     * @return quantileImpl1m(p) if p is in (0, 1)
     */
    [[nodiscard]] T Quantile1m(double p) const;

    /**
     * @param p
     * @param y
     * Fills vector y with Quantile(p)
     */
    void QuantileFunction(std::span<const double> p, std::span<T> y) const;

    /**
     * @param p
     * @param y
     * Fills vector y with Quantile1m(p)
     */
    void QuantileFunction1m(std::span<const double> p, std::span<T> y) const;

    /**
     * @param t
     * @return CFImpl for t != 0 and 1 otherwise
     */
    [[nodiscard]] std::complex<double> CF(double t) const;

    /**
     * @param t input vector
     * @param y output vector: y = CF(t)
     */
    void CharacteristicFunction(std::span<const double> t, std::span<std::complex<double>> y) const;

    /**
     * @param x input parameter
     * @return hazard function
     */
    [[nodiscard]] virtual double Hazard(const T &x) const = 0;

    /**
     * @param x input vector
     * @param y output vector: y = Hazard(x)
     */
    void HazardFunction(std::span<const T> x, std::span<double> y) const;

    /**
     * @return such x that F(x) = 0.5
     */
    [[nodiscard]] virtual T Median() const;

    /**
     * Mode is the value, which has the largest probability to happen,
     * or, in the case of continuous distribution, has the largest value
     * of density function.
     * @return mode. In case, when it's not unique, return any of them.
     */
    [[nodiscard]] virtual T Mode() const = 0;

    /**
     * @return E[((X - μ) / σ) ^ 3]
     * where μ is central moment and σ is standard deviation
     */
    [[nodiscard]] virtual long double Skewness() const;

    /**
     * @return unbiased kurtosis = μ_4 / σ ^ 4
     */
    virtual long double Kurtosis() const;

    /**
     * @return E[((X - μ) / σ) ^ 4]  - 3
     * (fourth moment around the mean divided by the square of the variance of the probability distribution minus 3)
     */
    [[nodiscard]] virtual long double ExcessKurtosis() const;

    /**
     * @return E[X^2]
     */
    [[nodiscard]] virtual long double SecondMoment() const;

    /**
     * @return E[X^3]
     */
    [[nodiscard]] virtual long double ThirdMoment() const;

    /**
     * @return E[X^4]
     */
    [[nodiscard]] virtual long double FourthMoment() const;

    /**
     * @param sample
     * @return likelihood function for given sample
     */
    [[nodiscard]] virtual double LikelihoodFunction(std::span<const T> sample) const = 0;

    /**
     * @param sample
     * @return logarithm of likelihood function for given sample
     */
    [[nodiscard]] virtual double LogLikelihoodFunction(std::span<const T> sample) const = 0;

protected:
    /**
     * @param value
     * @param sample
     * @return true if all elements in sample are not greater than given value
     */
    static bool allElementsAreNotGreaterThan(T value, std::span<const T> sample);

    /**
     * @param value
     * @param sample
     * @return true if all elements in sample are not smaller than given value
     */
    static bool allElementsAreNotSmallerThan(T value, std::span<const T> sample);

    /**
     * @param sample
     * @return true if all elements in sample are non-negative
     */
    static bool allElementsAreNonNegative(std::span<const T> sample);

    /**
     * @param sample
     * @return true if all elements in sample are positive
     */
    static bool allElementsArePositive(std::span<const T> sample);

public:
    /**
     * @param sample
     * @return sum of all elements in a sample
     */
    [[nodiscard]] static long double GetSampleSum(std::span<const T> sample);

    /**
     * @param sample
     * @return arithmetic average
     */
    [[nodiscard]] static long double GetSampleMean(std::span<const T> sample);

    /**
     * @param sample
     * @return arithmetic log-average
     */
    [[nodiscard]] static long double GetSampleLogMean(std::span<const T> sample);

    /**
     * @param sample
     * @param mean known mean value
     * @return sample second central moment
     */
    [[nodiscard]] static long double GetSampleVariance(std::span<const T> sample, double mean);

    /**
     * @param sample
     * @param logMean known log-mean value
     * @return sample log-variance
     */
    [[nodiscard]] static long double GetSampleLogVariance(std::span<const T> sample, double logMean);

    /**
     * @param sample
     * @return sample mean and variance
     */
    [[nodiscard]] static LongDoublePair GetSampleMeanAndVariance(std::span<const T> sample);

    /**
     * @param sample
     * @return sample log-mean and log-variance
     */
    [[nodiscard]] static LongDoublePair GetSampleLogMeanAndVariance(std::span<const T> sample);

    /**
     * @param sample
     * @return sample mean, variance, skewness and excess kurtosis
     */
    [[nodiscard]] static std::tuple<long double, long double, long double, long double>
    GetSampleStatistics(std::span<const T> sample);
};

} // namespace randlib

namespace randlib
{

template <typename T>
UnivariateDistribution<T>::UnivariateDistribution()
{
}

template <typename T>
bool UnivariateDistribution<T>::isLeftBounded() const
{
    Support supp = this->SupportType();
    return (supp == Support::RightSemiFinite || supp == Support::Finite);
}

template <typename T>
bool UnivariateDistribution<T>::isRightBounded() const
{
    Support supp = this->SupportType();
    return (supp == Support::LeftSemiFinite || supp == Support::Finite);
}

template <typename T>
T UnivariateDistribution<T>::Quantile(double p) const
{
    if (!(p >= 0.0 && p <= 1.0))
        throw std::invalid_argument("Probability p in quantile function should be in interval [0, 1]");
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

template <typename T>
T UnivariateDistribution<T>::Quantile1m(double p) const
{
    if (!(p >= 0.0 && p <= 1.0))
        throw std::invalid_argument("Probability p in quantile function should be in interval [0, 1]");
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

template <typename T>
void UnivariateDistribution<T>::QuantileFunction(std::span<const double> p, std::span<T> y) const
{
    const size_t size = std::min(p.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = this->Quantile(p[i]);
}

template <typename T>
void UnivariateDistribution<T>::QuantileFunction1m(std::span<const double> p, std::span<T> y) const
{
    const size_t size = std::min(p.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = this->Quantile1m(p[i]);
}

template <typename T>
std::complex<double> UnivariateDistribution<T>::CF(double t) const
{
    if (t == 0.0)
        return 1.0;
    return (t > 0.0) ? this->CFImpl(t) : std::conj(this->CFImpl(-t));
}

template <typename T>
std::complex<double> UnivariateDistribution<T>::CFImpl(double t) const
{
    T leftBound = this->MinValue(), rightBound = this->MaxValue();
    if (leftBound == rightBound)
        return std::complex<double>(std::cos(t * leftBound), std::sin(t * leftBound));

    double re = this->ExpectedValue([t](double x) { return std::cos(t * x); }, leftBound, rightBound);

    double im = this->ExpectedValue([t](double x) { return std::sin(t * x); }, leftBound, rightBound);

    return std::complex<double>(re, im);
}

template <typename T>
void UnivariateDistribution<T>::CharacteristicFunction(std::span<const double> t,
                                                       std::span<std::complex<double>> y) const
{
    const size_t size = std::min(t.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = this->CF(t[i]);
}

template <typename T>
void UnivariateDistribution<T>::HazardFunction(std::span<const T> x, std::span<double> y) const
{
    const size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = this->Hazard(x[i]);
}

template <typename T>
T UnivariateDistribution<T>::Median() const
{
    return quantileImpl(0.5);
}

template <typename T>
long double UnivariateDistribution<T>::Skewness() const
{
    long double var = this->Variance();
    if (!std::isfinite(var))
        return NAN;
    long double mu = this->Mean(); /// var is finite, so is mu

    long double sum = this->ExpectedValue(
        [mu](double x) {
            double xmmu = x - mu;
            double skewness = xmmu * xmmu * xmmu;
            return skewness;
        },
        this->MinValue(), this->MaxValue());

    return sum / std::pow(var, 1.5);
}

template <typename T>
long double UnivariateDistribution<T>::ExcessKurtosis() const
{
    long double var = this->Variance();
    if (!std::isfinite(var))
        return NAN;
    long double mu = this->Mean(); /// var is finite, so is mu

    long double sum = this->ExpectedValue(
        [mu](double x) {
            double xmmu = x - mu;
            double kurtosisSqrt = xmmu * xmmu;
            double kurtosis = kurtosisSqrt * kurtosisSqrt;
            return kurtosis;
        },
        this->MinValue(), this->MaxValue());

    return sum / (var * var) - 3;
}

template <typename T>
long double UnivariateDistribution<T>::Kurtosis() const
{
    return this->ExcessKurtosis() + 3.0;
}

template <typename T>
long double UnivariateDistribution<T>::SecondMoment() const
{
    long double mean = this->Mean();
    return mean * mean + this->Variance();
}

template <typename T>
long double UnivariateDistribution<T>::ThirdMoment() const
{
    long double mean = this->Mean();
    long double variance = this->Variance();
    long double skewness = this->Skewness();

    long double moment = skewness * std::sqrt(variance) * variance;
    moment += mean * mean * mean;
    moment += 3 * mean * variance;
    return moment;
}

template <typename T>
long double UnivariateDistribution<T>::FourthMoment() const
{
    long double mean = this->Mean();
    long double variance = this->Variance();
    long double moment3 = this->ThirdMoment();
    long double kurtosis = this->Kurtosis();
    long double meanSq = mean * mean;

    long double moment = kurtosis * variance * variance;
    moment -= 6 * meanSq * variance;
    moment -= 3 * meanSq * meanSq;
    moment += 4 * mean * moment3;
    return moment;
}

template <typename T>
bool UnivariateDistribution<T>::allElementsAreNotGreaterThan(T value, std::span<const T> sample)
{
    for (const T &var : sample) {
        if (var > value)
            return false;
    }
    return true;
}

template <typename T>
bool UnivariateDistribution<T>::allElementsAreNotSmallerThan(T value, std::span<const T> sample)
{
    for (const T &var : sample) {
        if (var < value)
            return false;
    }
    return true;
}

template <typename T>
bool UnivariateDistribution<T>::allElementsAreNonNegative(std::span<const T> sample)
{
    return allElementsAreNotSmallerThan(0, sample);
}

template <typename T>
bool UnivariateDistribution<T>::allElementsArePositive(std::span<const T> sample)
{
    for (const T &var : sample) {
        if (var <= 0)
            return false;
    }
    return true;
}

template <typename T>
long double UnivariateDistribution<T>::GetSampleSum(std::span<const T> sample)
{
    return std::accumulate(sample.begin(), sample.end(), 0.0);
}

template <typename T>
long double UnivariateDistribution<T>::GetSampleMean(std::span<const T> sample)
{
    size_t n = sample.size();
    return (n > 0) ? GetSampleSum(sample) / n : 0.0;
}

template <typename T>
long double UnivariateDistribution<T>::GetSampleLogMean(std::span<const T> sample)
{
    long double sum = 0.0;
    for (const T &var : sample)
        sum += std::log(var);
    return sum / sample.size();
}

template <typename T>
long double UnivariateDistribution<T>::GetSampleVariance(std::span<const T> sample, double mean)
{
    long double sum = 0.0l;
    for (const T &var : sample) {
        double temp = var - mean;
        sum += temp * temp;
    }
    return sum / sample.size();
}

template <typename T>
long double UnivariateDistribution<T>::GetSampleLogVariance(std::span<const T> sample, double logMean)
{
    long double sum = 0.0l;
    for (const T &var : sample) {
        double temp = std::log(var) - logMean;
        sum += temp * temp;
    }
    return sum / sample.size();
}

template <typename T>
LongDoublePair UnivariateDistribution<T>::GetSampleMeanAndVariance(std::span<const T> sample)
{
    /// Welford's stable method
    long double m = 0.0l, v = 0.0l;
    int n = sample.size();
    for (int i = 0; i < n; ++i) {
        double x = sample[i];
        double diff = x - m;
        m += diff / (i + 1);
        v += diff * (x - m);
    }
    return std::make_pair(m, v / n);
}

template <typename T>
LongDoublePair UnivariateDistribution<T>::GetSampleLogMeanAndVariance(std::span<const T> sample)
{
    /// Welford's stable method
    long double m = 0.0l, v = 0.0l;
    int n = sample.size();
    for (int i = 0; i < n; ++i) {
        double logX = std::log(sample[i]);
        double diff = logX - m;
        m += diff / (i + 1);
        v += diff * (logX - m);
    }
    return std::make_pair(m, v / n);
}

template <typename T>
std::tuple<long double, long double, long double, long double>
UnivariateDistribution<T>::GetSampleStatistics(std::span<const T> sample)
{
    /// Terriberry's extension for skewness and kurtosis
    long double M1{}, M2{}, M3{}, M4{};
    long double m1{}, m2{}, m3{}, m4{};
    size_t n = sample.size(), k = 0;
    size_t t = 0;
    static constexpr size_t BIG_NUMBER = 10000;
    for (const T &var : sample) {
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
            M2 += m2 + (DeltaSq * t * BIG_NUMBER) / tp1;
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
    long double res = static_cast<long double>(k) / BIG_NUMBER;
    if (res != 0) {
        long double Delta = m1 - M1;
        long double DeltaSq = Delta * Delta;
        long double tpres = t + res;
        /// M4
        M4 += m4 + (DeltaSq * DeltaSq * t * res * (t * t - res * t + res * res) * BIG_NUMBER) / (tpres * tpres * tpres);
        M4 += 6 * DeltaSq * (t * t * m2 + res * res * M2) / (tpres * tpres);
        M4 += 4 * Delta * (t * m3 - res * M3) / tpres;
        /// M3
        M3 += m3 + (DeltaSq * Delta * t * res * (t - res) * BIG_NUMBER) / (tpres * tpres);
        M3 += 3 * Delta * (t * m2 - res * M2) / tpres;
        /// M2 and M1
        M2 += m2 + (DeltaSq * t * res * BIG_NUMBER) / tpres;
        M1 += Delta * res / tpres;
    }

    long double variance = M2 / n;
    long double skewness = std::sqrt(n) * M3 / std::pow(M2, 1.5);
    long double exkurtosis = (n * M4) / (M2 * M2) - 3.0;
    return std::make_tuple(M1, variance, skewness, exkurtosis);
}

} // namespace randlib

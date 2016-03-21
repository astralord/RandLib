#include "NegativeBinomialRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

template< typename T >
NegativeBinomialRand<T>::NegativeBinomialRand(T number, double probability)
{
    setParameters(number, probability);
}

template< typename T >
std::string NegativeBinomialRand<T>::name()
{
    return "Negative Binomial(" + toStringWithPrecision(getNumber()) + ", " + toStringWithPrecision(getProbability()) + ")";
}

template< typename T >
void NegativeBinomialRand<T>::setValidParameters(T number, double probability)
{
    r = number;
    if (r <= 0.0)
        r = 1;

    p = probability;
    if (p > 1.0 || p < 0.0)
        p = 0.5;
    q = 1.0 - p;
}

template<>
void NegativeBinomialRand<int>::setParameters(int number, double probability)
{
    setValidParameters(number, probability);
    pdfCoef = std::pow(p, r) / RandMath::factorial(r - 1);

    if (r < 10)
    {
        /// we use two different generators for two different cases
        /// if p < 0.2 then the tail is too heavy
        /// (probability to be in main body is less than 0.977482)
        /// thus we return highest integer less than variate from exponential distribution
        /// otherwise we choose table method
        if (p < 0.2) {
            table[0] = -std::log(q);
        }
        else
        {
            table[0] = p;
            double prod = p;
            for (int i = 1; i < tableSize; ++i)
            {
                prod *= q;
                table[i] = table[i - 1] + prod;
            }
        }
    }
    else {
        Y.setParameters(r, q / p);
    }
}

template<>
void NegativeBinomialRand<double>::setParameters(double number, double probability)
{
    setValidParameters(number, probability);
    pdfCoef = std::pow(p, r) / std::tgamma(r);
    Y.setParameters(r, q / p);
}

template<>
double NegativeBinomialRand<double>::P(int k) const
{
    return (k < 0) ? 0 : pdfCoef * std::tgamma(r + k) / RandMath::factorial(k) * std::pow(q, k);
}

template<>
double NegativeBinomialRand<int>::P(int k) const
{
    return (k < 0) ? 0 : pdfCoef * RandMath::factorial(r + k - 1) / RandMath::factorial(k) * std::pow(q, k);
}

template< typename T >
double NegativeBinomialRand<T>::F(double x) const
{
    if (x < 0.0)
        return 0.0;
    return 1.0 - RandMath::regularizedBetaFun(q, std::floor(x) + 1, r);
}

template< typename T >
double NegativeBinomialRand<T>::variateThroughGammaPoisson() const
{
    return PoissonRand::variate(Y.variate());
}

template<>
double NegativeBinomialRand<double>::variate() const
{
    return variateThroughGammaPoisson();
}

template<>
double NegativeBinomialRand<int>::variateGeometricByTable() const
{
    double U = UniformRand::standardVariate();
    /// handle tail by recursion
    if (U > table[tableSize - 1])
        return tableSize + variateGeometricByTable();
    /// handle the main body
    int x = 0;
    while (U > table[x])
        ++x;
    return x;
}

template<>
double NegativeBinomialRand<int>::variateGeometricThroughExponential() const
{
    return std::floor(ExponentialRand::standardVariate() * table[0]);
}

template<>
double NegativeBinomialRand<int>::variateThroughGeometric() const
{
    double res = 0;
    if (p < 0.2)
    {
        for (int i = 0; i < r; ++i) {
            res += variateGeometricThroughExponential();
        }
    }
    else
    {
        for (int i = 0; i < r; ++i) {
            res += variateGeometricByTable();
        }
    }
    return res;
}

template<>
double NegativeBinomialRand<int>::variate() const
{
    if (r < 10)
        return variateThroughGeometric();
    return variateThroughGammaPoisson();
}

template<>
void NegativeBinomialRand<double>::sample(QVector<double> &outputData) const
{
    for (double &var : outputData)
        var = variateThroughGammaPoisson();
}

template<>
void NegativeBinomialRand<int>::sample(QVector<double> &outputData) const
{
    if (r < 10)
    {
        if (p < 0.2)
        {
            for (double &var : outputData)
            {
                var = 0;
                for (int i = 0; i < r; ++i)
                    var += variateGeometricThroughExponential();
            }
        }
        else
        {
            for (double &var : outputData)
            {
                var = 0;
                for (int i = 0; i < r; ++i)
                    var += variateGeometricByTable();
            }
        }
    }
    else
    {
        for (double &var : outputData)
            var = variateThroughGammaPoisson();
    }
}

template< typename T >
double NegativeBinomialRand<T>::Mean() const
{
    return q * r / p;
}

template< typename T >
double NegativeBinomialRand<T>::Variance() const
{
    return q * r / (p * p);
}

template< typename T >
std::complex<double> NegativeBinomialRand<T>::CF(double t) const
{
    std::complex<double> it(0, t);
    std::complex<double> denominator = 1.0 - q * std::exp(it);
    return std::pow(p / denominator, r);
}

template< typename T >
double NegativeBinomialRand<T>::Mode() const
{
    return (r > 1) ? std::floor((r - 1) * q / p) : 0;
}

template< typename T >
double NegativeBinomialRand<T>::Skewness() const
{
    return (1 + q) / std::sqrt(q * r);
}

template< typename T >
double NegativeBinomialRand<T>::ExcessKurtosis() const
{
    double kurtosis = p * p;
    kurtosis /= q;
    kurtosis += 6;
    return kurtosis / r;
}


template class NegativeBinomialRand<int>;
template class NegativeBinomialRand<double>;

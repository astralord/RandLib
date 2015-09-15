#include "NegativeBinomialRand.h"

template< typename T >
NegativeBinomialRand<T>::NegativeBinomialRand(T number, double probability) :
    G(probability)
{
    setParameters(number, probability);
}

template< typename T >
std::string NegativeBinomialRand<T>::name()
{
    return "Negative Binomial(" + toStringWithPrecision(getNumber()) + ", " + toStringWithPrecision(getProbability()) + ")";
}

template< typename T >
void NegativeBinomialRand<T>::setParameters(T number, double probability)
{
    r = std::max(number, static_cast<T>(1.0));

    p = std::min(std::max(probability, MIN_POSITIVE), 1.0);
    G.setProbability(1 - p);

    pdfCoef = std::pow(1 - p, r) / std::tgamma(r);
    Y.setParameters(r, p / (1 - p));
}

template <>
double NegativeBinomialRand<double>::P(int k) const
{
    return (k < 0) ? 0 : pdfCoef * std::tgamma(r + k) / RandMath::factorial(k) * std::pow(p, k);
}

template <>
double NegativeBinomialRand<int>::P(int k) const
{
    return (k < 0) ? 0 : pdfCoef * RandMath::factorial(r + k - 1) / RandMath::factorial(k) * std::pow(p, k);
}

template< typename T >
double NegativeBinomialRand<T>::F(double x) const
{
    return 1.0 - RandMath::incompleteBetaFun(p, std::floor(x) + 1, r);
}

template<>
double NegativeBinomialRand<double>::variate() const
{
    return variateThroughGammaPoisson();
}

template<>
double NegativeBinomialRand<int>::variate() const
{
    if (r < 10)
        return variateThroughGeometric();
    return variateThroughGammaPoisson();
}

template< typename T >
double NegativeBinomialRand<T>::variateThroughGeometric() const
{
    double res = 0;
    for (int i = 0; i != static_cast<int>(r); ++i)
        res += G.variate();
    return res;
}

template< typename T >
double NegativeBinomialRand<T>::variateThroughGammaPoisson() const
{
    return PoissonRand::variate(Y.variate());
}

template< typename T >
double NegativeBinomialRand<T>::Mean() const
{
    return p * r / (1 - p);
}

template< typename T >
double NegativeBinomialRand<T>::Variance() const
{
    return Mean() / (1 - p);
}

template< typename T >
std::complex<double> NegativeBinomialRand<T>::CF(double t) const
{
    double numerator = 1 - p;
    std::complex denominator(0, t);
    denominator = 1 - p * std::exp(denominator);
    return std::pow(numerator / denominator, r);
}

template< typename T >
double NegativeBinomialRand<T>::Mode() const
{
    return (r > 1) ? std::floor((r - 1) * p / (1 - p)) : 0;
}

template< typename T >
double NegativeBinomialRand<T>::Skewness() const
{
    return (1 + p) / std::sqrt(p * r);
}

template< typename T >
double NegativeBinomialRand<T>::ExcessKurtosis() const
{
    double kurtosis = (1 - p);
    kurtosis *= kurtosis;
    kurtosis /= p;
    kurtosis += 6;
    return kurtosis / r;
}

template class NegativeBinomialRand<int>;
template class NegativeBinomialRand<double>;

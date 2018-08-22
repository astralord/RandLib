#include "UniformDiscreteRand.h"

template< typename IntType >
UniformDiscreteRand<IntType>::UniformDiscreteRand(IntType minValue, IntType maxValue)
{
    SetBoundaries(minValue, maxValue);
}

template< typename IntType >
String UniformDiscreteRand<IntType>::Name() const
{
    return "Uniform Discrete(" + this->toStringWithPrecision(MinValue()) + ", " + this->toStringWithPrecision(MaxValue()) + ")";
}

template< typename IntType >
void UniformDiscreteRand<IntType>::SetBoundaries(IntType minValue, IntType maxValue)
{
    if (minValue > maxValue)
        throw std::invalid_argument("Uniform discrete distribution: minimal value shouldn't be larger than maximum value");

    a = minValue;
    b = maxValue;

    n = b - a + 1;
    nInv = 1.0 / n;
    logN = std::log(n);

    unsigned long long MAX_RAND = this->localRandGenerator.MaxValue();
    MAX_RAND_UNBIASED = MAX_RAND - MAX_RAND % n - 1;
}

template< typename IntType >
double UniformDiscreteRand<IntType>::P(const IntType & k) const
{
    return (k < a || k > b) ? 0.0 : nInv;
}

template< typename IntType >
double UniformDiscreteRand<IntType>::logP(const IntType & k) const
{
    return (k < a || k > b) ? -INFINITY : -logN;
}

template< typename IntType >
double UniformDiscreteRand<IntType>::F(const IntType & k) const
{
    if (k < a)
        return 0.0;
    if (k > b)
        return 1.0;
    return (k - a + 1) * nInv;
}

template< typename IntType >
IntType UniformDiscreteRand<IntType>::Variate() const
{
    unsigned long intVar;
    do {
        intVar = this->localRandGenerator.Variate();
    } while (intVar > MAX_RAND_UNBIASED);
    return a + (intVar % n);
}

template< typename IntType >
IntType UniformDiscreteRand<IntType>::StandardVariate(IntType minValue, IntType maxValue, RandGenerator &randGenerator)
{
    unsigned long long MAX_RAND = randGenerator.MaxValue();
    IntType n = maxValue - minValue + 1;
    if (n <= 1)
        return minValue;
    unsigned long long MAX_RAND_UNBIASED = MAX_RAND - MAX_RAND % n - 1;
    unsigned long intVar;
    do {
        intVar = randGenerator.Variate();
    } while (intVar > MAX_RAND_UNBIASED);
    return minValue + (intVar % n);
}

template< typename IntType >
long double UniformDiscreteRand<IntType>::Mean() const
{
    return 0.5 * (b + a);
}

template< typename IntType >
long double UniformDiscreteRand<IntType>::Variance() const
{
    double nm1 = n - 1;
    double np1 = n + 1;
    return nm1 * np1 / 12;
}

template< typename IntType >
IntType UniformDiscreteRand<IntType>::Median() const
{
    return (b + a) >> 1;
}

template< typename IntType >
IntType UniformDiscreteRand<IntType>::Mode() const
{
    /// this can be any value in [a, b]
    return 0.5 * (a + b);
}

template< typename IntType >
long double UniformDiscreteRand<IntType>::Skewness() const
{
    return 0.0;
}

template< typename IntType >
long double UniformDiscreteRand<IntType>::ExcessKurtosis() const
{
    double kurt = n;
    kurt *= n;
    --kurt;
    kurt = 2.0 / kurt;
    ++kurt;
    return -1.2 * kurt;
}

template< typename IntType >
std::complex<double> UniformDiscreteRand<IntType>::CFImpl(double t) const
{
    double at = a * t;
    double bp1t = (b + 1) * t;
    double reNum = std::cos(at) - std::cos(bp1t);
    double imNum = std::sin(at) - std::sin(bp1t);
    std::complex<double> numerator(reNum, imNum);
    std::complex<double> denominator(1.0 - std::cos(t), -std::sin(t));
    return nInv * numerator / denominator;
}

template< typename IntType >
double UniformDiscreteRand<IntType>::Entropy() const
{
    return logN;
}

template< typename IntType >
double UniformDiscreteRand<IntType>::LikelihoodFunction(const std::vector<IntType> &sample) const
{
    bool sampleIsInsideInterval = this->allElementsAreNotSmallerThan(a, sample) && this->allElementsAreNotLargerThan(b, sample);
    return sampleIsInsideInterval ? std::pow(n, -sample.size()) : 0.0;
}

template< typename IntType >
double UniformDiscreteRand<IntType>::LogLikelihoodFunction(const std::vector<IntType> &sample) const
{
    bool sampleIsInsideInterval = this->allElementsAreNotSmallerThan(a, sample) && this->allElementsAreNotLargerThan(b, sample);
    int sample_size = sample.size();
    return sampleIsInsideInterval ? -sample_size * logN : -INFINITY;
}

template< typename IntType >
void UniformDiscreteRand<IntType>::Fit(const std::vector<IntType> &sample)
{
    IntType minVar = *std::min_element(sample.begin(), sample.end());
    IntType maxVar = *std::max_element(sample.begin(), sample.end());
    this->SetBoundaries(minVar, maxVar);
}


template class UniformDiscreteRand<int>;
template class UniformDiscreteRand<long int>;
template class UniformDiscreteRand<long long int>;

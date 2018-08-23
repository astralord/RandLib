#include "HyperGeometricRand.h"

template < typename IntType >
HyperGeometricRand<IntType>::HyperGeometricRand(IntType totalSize, IntType drawsNum, IntType successesNum)
{
    SetParameters(totalSize, drawsNum, successesNum);
}

template < typename IntType >
String HyperGeometricRand<IntType>::Name() const
{
    return "Hypergeometric(" + this->toStringWithPrecision(N) + ", "
                             + this->toStringWithPrecision(n) + ", "
                             + this->toStringWithPrecision(K) + ")";
}

template < typename IntType >
void HyperGeometricRand<IntType>::SetParameters(IntType totalSize, IntType drawsNum, IntType successesNum)
{
    if (totalSize <= 0 || drawsNum <= 0 || successesNum <= 0)
        throw std::invalid_argument("HyperGeometric distribution: all parameters should be positive");
    if (drawsNum > totalSize)
        throw std::invalid_argument("HyperGeometric distribution: total size should be greater than draws number");
    if (successesNum > totalSize)
        throw std::invalid_argument("HyperGeometric distribution: total size should be greater than successes number");

    N = totalSize;
    n = drawsNum;
    K = successesNum;

    p0 = static_cast<double>(K) / N;
    pmfCoef = RandMath::lfact(K);
    pmfCoef += RandMath::lfact(N - K);
    pmfCoef += RandMath::lfact(N - n);
    pmfCoef += RandMath::lfact(n);
    pmfCoef -= RandMath::lfact(N);
}

template < typename IntType >
double HyperGeometricRand<IntType>::P(const IntType & k) const
{
    return (k < MinValue() || k > MaxValue()) ? 0.0 : std::exp(logP(k));
}

template < typename IntType >
double HyperGeometricRand<IntType>::logP(const IntType & k) const
{
    if (k < MinValue() || k > MaxValue())
        return -INFINITY;
    double y = RandMath::lfact(k);
    y += RandMath::lfact(K - k);
    y += RandMath::lfact(n - k);
    y += RandMath::lfact(N - K - n + k);
    return pmfCoef - y;
}

template < typename IntType >
double HyperGeometricRand<IntType>::F(const IntType & k) const
{
    if (k < MinValue())
        return 0.0;
    IntType maxVal = MaxValue();
    if (k >= maxVal)
        return 1.0;
    if (k <= 0.5 * maxVal) {
        /// sum P(X = i) going forward until k
        double sum = 0;
        for (IntType i = 0; i <= k; ++i)
            sum += P(i);
        return sum;
    }
    /// going backwards is faster
    double sum = 1.0;
    for (IntType i = k + 1; i <= maxVal; ++i)
        sum -= P(i);
    return sum;
}

template < typename IntType >
IntType HyperGeometricRand<IntType>::Variate() const
{
    double p = p0;
    IntType sum = 0;
    for (int i = 1; i <= n; ++i)
    {
        if (BernoulliRand::Variate(p, this->localRandGenerator) && ++sum >= K)
            return sum;
        p = K - sum;
        p /= N - i;
    }
    return sum;
}

template < typename IntType >
long double HyperGeometricRand<IntType>::Mean() const
{
    return static_cast<double>(n * K) / N;
}

template < typename IntType >
long double HyperGeometricRand<IntType>::Variance() const
{
    long double numerator = n;
    numerator *= K;
    numerator *= N - K;
    numerator *= N  - n;
    long double denominator = N;
    denominator *= N;
    denominator *= N - 1;
    return numerator / denominator;
}

template < typename IntType >
IntType HyperGeometricRand<IntType>::Mode() const
{
    double mode = (n + 1) * (K + 1);
    return std::floor(mode / (N + 2));
}

template < typename IntType >
long double HyperGeometricRand<IntType>::Skewness() const
{
    long double skewness = N - 1;
    skewness /= n;
    skewness /= K;
    skewness /= N - K;
    skewness /= N - n;
    skewness = std::sqrt(skewness);
    skewness *= N - 2 * K;
    skewness *= N - 2 * n;
    return skewness / (N - 2);
}

template < typename IntType >
long double HyperGeometricRand<IntType>::ExcessKurtosis() const
{
    long double numerator = N;
    numerator *= (N + 1);
    long double a1 = K;
    a1 *= N - K;
    long double a2 = n;
    a2 *= N - n;
    numerator -= 6 * (a1 + a2);
    numerator *= N;
    numerator *= N;
    numerator *= N - 1;
    long double aux = n;
    aux *= K;
    aux *= N - K;
    aux *= N - n;
    numerator += 6 * aux * (5 * N - 6);
    long double denominator = aux * (N - 2) * (N - 3);
    return numerator / denominator;
}

template class HyperGeometricRand<int>;
template class HyperGeometricRand<long int>;
template class HyperGeometricRand<long long int>;

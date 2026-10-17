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
    if (!(totalSize > 0 && drawsNum > 0 && successesNum > 0))
        throw std::invalid_argument("HyperGeometric distribution: all parameters should be positive");
    if (!(drawsNum <= totalSize))
        throw std::invalid_argument("HyperGeometric distribution: total size should be greater than draws number");
    if (!(successesNum <= totalSize))
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
    IntType minVal = MinValue();
    if (k < minVal)
        return 0.0;
    IntType maxVal = MaxValue();
    if (k >= maxVal)
        return 1.0;
    if (k - minVal <= maxVal - k) {
        /// sum P(X = i) going forward until k
        double sum = 0;
        for (IntType i = minVal; i <= k; ++i)
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
    for (IntType i = 1; i <= n; ++i)
    {
        if (BernoulliRand::Variate(p, this->localRandGenerator) && ++sum >= K)
            return sum;
        if (i < n) {
            p = K - sum;
            p /= N - i;
        }
    }
    return sum;
}

template < typename IntType >
long double HyperGeometricRand<IntType>::Mean() const
{
    return static_cast<long double>(n) * K / N;
}

template < typename IntType >
long double HyperGeometricRand<IntType>::Variance() const
{
    if (n == N || K == N)
        return 0.0L;
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
    long double mode = (static_cast<long double>(n) + 1.0L) * (static_cast<long double>(K) + 1.0L);
    return static_cast<IntType>(std::floor(mode / (static_cast<long double>(N) + 2.0L)));
}

template < typename IntType >
long double HyperGeometricRand<IntType>::Skewness() const
{
    if (n == N || K == N)
        return NAN;
    if (N == 2)
        return 0.0L;
    long double population = N, draws = n, successes = K;
    long double skewness = population - 1.0L;
    skewness /= draws;
    skewness /= successes;
    skewness /= population - successes;
    skewness /= population - draws;
    skewness = std::sqrt(skewness);
    skewness *= population - 2.0L * successes;
    skewness *= population - 2.0L * draws;
    return skewness / (population - 2.0L);
}

template < typename IntType >
long double HyperGeometricRand<IntType>::ExcessKurtosis() const
{
    long double variance = Variance();
    if (variance == 0.0L)
        return NAN;
    if (N <= 3) {
        long double mean = Mean(), fourthMoment = 0.0L;
        for (IntType k = MinValue(); k <= MaxValue(); ++k) {
            long double centered = k - mean;
            fourthMoment += centered * centered * centered * centered * P(k);
        }
        return fourthMoment / (variance * variance) - 3.0L;
    }

    long double population = N, draws = n, successes = K;
    long double numerator = population * (population + 1.0L);
    long double a1 = successes * (population - successes);
    long double a2 = draws * (population - draws);
    numerator -= 6 * (a1 + a2);
    numerator *= population * population * (population - 1.0L);
    long double aux = draws * successes * (population - successes) * (population - draws);
    numerator += 6.0L * aux * (5.0L * population - 6.0L);
    long double denominator = aux * (population - 2.0L) * (population - 3.0L);
    return numerator / denominator;
}

template class HyperGeometricRand<int>;
template class HyperGeometricRand<long int>;
template class HyperGeometricRand<long long int>;

#include "NegativeHyperGeometricRand.h"

template < typename IntType >
NegativeHyperGeometricRand<IntType>::NegativeHyperGeometricRand(IntType totalSize, IntType totalSuccessesNum, IntType limitSuccessesNum)
{
    SetParameters(totalSize, totalSuccessesNum, limitSuccessesNum);
}

template < typename IntType >
String NegativeHyperGeometricRand<IntType>::Name() const
{
    return "Negative hypergeometric(" + this->toStringWithPrecision(N) + ", "
                                      + this->toStringWithPrecision(M) + ", "
                                      + this->toStringWithPrecision(m) + ")";
}

template < typename IntType >
void NegativeHyperGeometricRand<IntType>::SetParameters(IntType totalSize, IntType totalSuccessesNum, IntType limitSuccessesNum)
{
    if (totalSize <= 0 || totalSuccessesNum <= 0 || limitSuccessesNum <= 0)
        throw std::invalid_argument("Negative-HyperGeometric distribution: all parameters should be positive");
    if (totalSuccessesNum > totalSize)
        throw std::invalid_argument("Negative-HyperGeometric distribution: total size shouldn't be smaller than total successes number");
    if (limitSuccessesNum > totalSuccessesNum)
        throw std::invalid_argument("Negative-HyperGeometric distribution: total successes number shouldn't be smaller than limit successes number");

    N = totalSize;
    M = totalSuccessesNum;
    m = limitSuccessesNum;

    p0 = static_cast<double>(M) / N;
    pmfCoef = RandMath::lfact(M);
    pmfCoef += RandMath::lfact(N - M);
    pmfCoef -= RandMath::lfact(m - 1);
    pmfCoef -= RandMath::lfact(M - m);
    pmfCoef -= RandMath::lfact(N);
}

template < typename IntType >
double NegativeHyperGeometricRand<IntType>::P(const IntType & k) const
{
    return (k < MinValue() || k > MaxValue()) ? 0.0 : std::exp(logP(k));
}

template < typename IntType >
double NegativeHyperGeometricRand<IntType>::logP(const IntType & k) const
{
    if (k < MinValue() || k > MaxValue())
        return -INFINITY;
    double p = RandMath::lfact(k + m - 1);
    p += RandMath::lfact(N - m - k);
    p -= RandMath::lfact(k);
    p -= RandMath::lfact(N - M - k);
    return p + pmfCoef;
}

template < typename IntType >
double NegativeHyperGeometricRand<IntType>::F(const IntType &k) const
{
    // relation with hypergeometric distribution can be used here instead
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
IntType NegativeHyperGeometricRand<IntType>::Variate() const
{
    double p = p0;
    IntType successesNum = 0;
    IntType num = 0;
    while (successesNum < m) {
        ++num;
        if (BernoulliRand::Variate(p, this->localRandGenerator) && ++successesNum == num - N + M)
            return N - M;
        p = M - successesNum;
        p /= N - num;
    }
    return num - successesNum;
}

template < typename IntType >
long double NegativeHyperGeometricRand<IntType>::Mean() const
{
    double mean = m;
    mean *= N - M;
    mean /= M + 1;
    return mean;
}

template < typename IntType >
long double NegativeHyperGeometricRand<IntType>::Variance() const
{
    double Mp1 = M + 1;
    double var = 1 - m / Mp1;
    var *= N + 1;
    var *= N - M;
    var /= Mp1 * (Mp1 + 1);
    return m * var;
}

template class NegativeHyperGeometricRand<int>;
template class NegativeHyperGeometricRand<long int>;
template class NegativeHyperGeometricRand<long long int>;

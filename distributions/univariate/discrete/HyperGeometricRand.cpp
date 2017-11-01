#include "HyperGeometricRand.h"

HyperGeometricRand::HyperGeometricRand(int totalSize, int drawsNum, int successesNum)
{
    SetParameters(totalSize, drawsNum, successesNum);
}

String HyperGeometricRand::Name() const
{
    return "Hypergeometric(" + toStringWithPrecision(N) + ", "
                             + toStringWithPrecision(n) + ", "
                             + toStringWithPrecision(K) + ")";
}

void HyperGeometricRand::SetParameters(int totalSize, int drawsNum, int successesNum)
{
    if (totalSize <= 0 || drawsNum <= 0 || successesNum <= 0)
        throw std::invalid_argument("HyperGeometric distribution: all parameters should be positive");
    if (drawsNum > totalSize)
        throw std::invalid_argument("HyperGeometric distribution: total size should be larger than draws number");
    if (successesNum > totalSize)
        throw std::invalid_argument("HyperGeometric distribution: total size should be larger than successes number");

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

double HyperGeometricRand::P(const int & k) const
{
    return (k < MinValue() || k > MaxValue()) ? 0.0 : std::exp(logP(k));
}

double HyperGeometricRand::logP(const int & k) const
{
    if (k < MinValue() || k > MaxValue())
        return -INFINITY;
    double y = RandMath::lfact(k);
    y += RandMath::lfact(K - k);
    y += RandMath::lfact(n - k);
    y += RandMath::lfact(N - K - n + k);
    return pmfCoef - y;
}

double HyperGeometricRand::F(const int & k) const
{
    if (k < MinValue())
        return 0.0;
    int maxVal = MaxValue();
    if (k >= maxVal)
        return 1.0;
    if (k <= 0.5 * maxVal) {
        /// sum P(X = i) going forward until k
        double sum = 0;
        for (int i = 0; i <= k; ++i)
            sum += P(i);
        return sum;
    }
    /// going backwards is faster
    double sum = 1.0;
    for (int i = k + 1; i <= maxVal; ++i)
        sum -= P(i);
    return sum;
}

int HyperGeometricRand::Variate() const
{
    double p = p0;
    int sum = 0;
    for (int i = 1; i <= n; ++i)
    {
        if (BernoulliRand::Variate(p) && ++sum >= K)
            return sum;
        p = K - sum;
        p /= N - i;
    }
    return sum;
}

double HyperGeometricRand::Mean() const
{
    return static_cast<double>(n * K) / N;
}

double HyperGeometricRand::Variance() const
{
    double numerator = n;
    numerator *= K;
    numerator *= N - K;
    numerator *= N  - n;
    double denominator = N;
    denominator *= N;
    denominator *= N - 1;
    return numerator / denominator;
}

int HyperGeometricRand::Mode() const
{
    double mode = (n + 1) * (K + 1);
    return std::floor(mode / (N + 2));
}

double HyperGeometricRand::Skewness() const
{
    double skewness = N - 1;
    skewness /= n;
    skewness /= K;
    skewness /= N - K;
    skewness *= N - n;
    skewness = std::sqrt(skewness);
    skewness *= N - 2 * K;
    skewness *= N - 2 * n;
    return skewness / (N - 2);
}

double HyperGeometricRand::ExcessKurtosis() const
{
    double numerator = N;
    numerator *= (N + 1);
    double a1 = K;
    a1 *= N - K;
    double a2 = n;
    a2 *= N - n;
    numerator -= 6 * (a1 + a2);
    numerator *= N;
    numerator *= N;
    numerator *= N - 1;
    double aux = n;
    aux *= K;
    aux *= N - K;
    aux *= N - n;
    numerator += 6 * aux * (5 * N - 6);
    double denominator = aux * (N - 2) * (N - 3);
    return numerator / denominator;
}


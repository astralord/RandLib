#include "HyperGeometricRand.h"

HyperGeometricRand::HyperGeometricRand(int totalSize, int drawsNum, int successesNum)
{
    SetParameters(totalSize, drawsNum, successesNum);
}

std::string HyperGeometricRand::Name() const
{
    return "Hypergeometric(" + toStringWithPrecision(N) + ", "
                             + toStringWithPrecision(n) + ", "
                             + toStringWithPrecision(K) + ")";
}

void HyperGeometricRand::SetParameters(int totalSize, int drawsNum, int successesNum)
{
    N = std::max(totalSize, 0);

    n = std::max(drawsNum, 0);
    n = std::min(N, n);

    K = std::max(successesNum, 0);
    K = std::min(N, K);

    p0 = static_cast<double>(K) / N;
    pmfDenominator = 1.0 / RandMath::binomialCoef(N, n);
}

double HyperGeometricRand::P(int k) const
{
    if (k < MinValue() || k > MaxValue())
        return 0.0;
    return RandMath::binomialCoef(K, k) * RandMath::binomialCoef(N - K, n - k) * pmfDenominator;
}

double HyperGeometricRand::F(int k) const
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
    double numerator = n * K * (N - K) * (N  - n);
    double denominator = N * N * (N - 1);
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
    skewness /= n * K * (N - K) * (N - n);
    skewness = std::sqrt(skewness);
    skewness *= N - K - K;
    skewness *= N - n - n;
    return skewness / (N - 2);
}

double HyperGeometricRand::ExcessKurtosis() const
{
    double numerator = N * (N + 1);
    numerator -= 6 * K * (N - K);
    numerator -= 6 * n * (N - n);
    numerator *= N * N * (N - 1);
    double aux = n * K * (N - K) * (N - n);
    numerator += 6 * aux * (5 * N - 6);
    double denominator = aux * (N - 2) * (N - 3);
    return numerator / denominator;
}


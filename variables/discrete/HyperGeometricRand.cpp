#include "HyperGeometricRand.h"

HyperGeometricRand::HyperGeometricRand(int totalSize, int drawsNum, int successesNum)
{
    setParameters(totalSize, drawsNum, successesNum);
}

std::string HyperGeometricRand::name()
{
    return "HyperGeometric(" + toStringWithPrecision(N) + ", "
                             + toStringWithPrecision(n) + ", "
                             + toStringWithPrecision(K) + ")";
}

void HyperGeometricRand::setParameters(int totalSize, int drawsNum, int successesNum)
{
    N = std::max(totalSize, 1);

    n = std::max(drawsNum, 1);
    n = std::min(N, n);

    K = std::max(successesNum, 1);
    K = std::min(N, K);

    p0 = static_cast<double>(K) / N;
    pdfDenominator = 1.0 / RandMath::binomialCoef(N, n);
}

double HyperGeometricRand::P(int k) const
{
    if (k < 0 || k > n || k > K || n - k > N - K)
        return 0;
    return RandMath::binomialCoef(K, k) * RandMath::binomialCoef(N - K, n - k) * pdfDenominator;
}

double HyperGeometricRand::F(double x) const
{
    double k = std::floor(x);
    if (k < 0 || n - k > N - K)
        return 0.0;
    if (k > n || k > K)
        return 1.0;
    double sum = 0;
    for (int i = 0; i <= k; ++i)
        sum += P(i);
    return sum;
}

double HyperGeometricRand::variate() const
{
    double p = p0, sum = 0;
    for (int i = 1; i <= n; ++i)
    {
        if (BernoulliRand::variate(p) && ++sum >= K)
            return sum;
        p = (K - sum) / (N - i);
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

double HyperGeometricRand::Mode() const
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


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
    N = totalSize;
    n = std::min(N, drawsNum);
    K = std::min(N, successesNum);

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
    int k = static_cast<int>(std::floor(x));
    double sum = 0;
    for (int i = 0; i < k; ++i)
        sum += P(i);
    return sum;
}

double HyperGeometricRand::variate() const
{
    size_t sum = 0;
    double p = static_cast<double>(K) / N;
    double successesLeft = K;
    for (int i = 0; i != n; ++i)
    {
        bool isSuccess = BernoulliRand::variate(p);
        sum += isSuccess;
        if (isSuccess)
            --successesLeft;
        if (successesLeft == 0)
            return sum;
        p = successesLeft / (N - i - 1);
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


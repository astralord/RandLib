#include "HyperGeometricRand.h"

HyperGeometricRand::HyperGeometricRand(size_t totalSize, size_t drawsNum, size_t successesNum)
{
    setParameters(totalSize, drawsNum, successesNum);
}

std::string HyperGeometricRand::name()
{
    return "HyperGeometric(" + toStringWithPrecision(N) + ", "
                             + toStringWithPrecision(n) + ", "
                             + toStringWithPrecision(K) + ")";
}

void HyperGeometricRand::setParameters(size_t totalSize, size_t drawsNum, size_t successesNum)
{
    N = totalSize;
    n = std::min(N, drawsNum);
    K = std::min(N, successesNum);

    pdfDenominator = 1.0 / RandMath::binomialCoef(N, n);
}

double HyperGeometricRand::P(int k) const
{
    size_t kUns = static_cast<size_t>(k);
    if (k < 0 || kUns > n || kUns > K || n - kUns > N - K)
        return 0;
    return RandMath::binomialCoef(K, kUns) * RandMath::binomialCoef(N - K, n - kUns) * pdfDenominator;
}

double HyperGeometricRand::F(double x) const
{
    size_t k = static_cast<size_t>(std::floor(x));
    double sum = 0;
    for (size_t i = 0; i < k; ++i)
        sum += P(i);
    return sum;
}

double HyperGeometricRand::variate() const
{
    size_t sum = 0;
    double p = static_cast<double>(K) / N;
    double successesLeft = K;
    for (size_t i = 0; i != n; ++i)
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

double HyperGeometricRand::E() const
{
    return static_cast<double>(n * K) / N;
}

double HyperGeometricRand::Var() const
{
    size_t numerator = n * K * (N - K) * (N  - n);
    size_t denominator = N * N * (N - 1);
    return static_cast<double>(numerator) / denominator;

}


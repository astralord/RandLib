#include "NegativeHyperGeometricRand.h"

NegativeHyperGeometricRand::NegativeHyperGeometricRand(int totalSize, int totalSuccessesNum, int limitSuccessesNum)
{
    SetParameters(totalSize, totalSuccessesNum, limitSuccessesNum);
}

String NegativeHyperGeometricRand::Name() const
{
    return "Negative hypergeometric(" + toStringWithPrecision(N) + ", "
                                      + toStringWithPrecision(M) + ", "
                                      + toStringWithPrecision(m) + ")";
}

void NegativeHyperGeometricRand::SetParameters(int totalSize, int totalSuccessesNum, int limitSuccessesNum)
{
    if (totalSize <= 0 || totalSuccessesNum <= 0 || limitSuccessesNum <= 0)
        throw std::invalid_argument("Negative-HyperGeometric distribution: all parameters should be positive");
    if (totalSuccessesNum > totalSize)
        throw std::invalid_argument("Negative-HyperGeometric distribution: total size should be larger than total successes number");
    if (limitSuccessesNum > totalSuccessesNum)
        throw std::invalid_argument("Negative-HyperGeometric distribution: total successes number should be larger than limit successes number");

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

double NegativeHyperGeometricRand::P(const int & k) const
{
    return (k < MinValue() || k > MaxValue()) ? 0.0 : std::exp(logP(k));
}

double NegativeHyperGeometricRand::logP(const int & k) const
{
    if (k < MinValue() || k > MaxValue())
        return -INFINITY;
    double p = RandMath::lfact(k + m - 1);
    p += RandMath::lfact(N - m - k);
    p -= RandMath::lfact(k);
    p -= RandMath::lfact(N - M - k);
    return p + pmfCoef;
}

double NegativeHyperGeometricRand::F(const int & k) const
{
    // relation with hypergeometric distribution can be used here instead
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

int NegativeHyperGeometricRand::Variate() const
{
    double p = p0;
    int successesNum = 0;
    int num = 0;
    while (successesNum < m) {
        ++num;
        if (BernoulliRand::Variate(p) && ++successesNum == num - N + M)
            return N - M;
        p = M - successesNum;
        p /= N - num;
    }
    return num - successesNum;
}

double NegativeHyperGeometricRand::Mean() const
{
    double mean = m;
    mean *= N - M;
    mean /= M + 1;
    return mean;
}

double NegativeHyperGeometricRand::Variance() const
{
    double Mp1 = M + 1;
    double var = 1 - m / Mp1;
    var *= N + 1;
    var *= N - M;
    var /= Mp1 * (Mp1 + 1);
    return m * var;
}

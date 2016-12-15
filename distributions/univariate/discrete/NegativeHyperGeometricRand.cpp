#include "NegativeHyperGeometricRand.h"

NegativeHyperGeometricRand::NegativeHyperGeometricRand(int totalSize, int totalSuccessesNum, int limitSuccessesNum)
{
    SetParameters(totalSize, totalSuccessesNum, limitSuccessesNum);
}

std::string NegativeHyperGeometricRand::Name() const
{
    return "Negative hypergeometric(" + toStringWithPrecision(N) + ", "
                                      + toStringWithPrecision(M) + ", "
                                      + toStringWithPrecision(m) + ")";
}

void NegativeHyperGeometricRand::SetParameters(int totalSize, int totalSuccessesNum, int limitSuccessesNum)
{
    N = std::max(0, totalSize);

    M = std::max(0, totalSuccessesNum);
    M = std::min(N, M);

    m = std::max(0, limitSuccessesNum);
    m = std::min(M, m);

    p0 = static_cast<double>(M) / N;
    pmfDenominator = 1.0 / RandMath::binomialCoef(N, M);
}

double NegativeHyperGeometricRand::P(int k) const
{
    if (k < MinValue() || k > MaxValue())
        return 0.0;
    double p = RandMath::binomialCoef(k + m - 1, k);
    p *= RandMath::binomialCoef(N - m - k, M - m);
    return p * pmfDenominator;
}

double NegativeHyperGeometricRand::F(int k) const
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
    return m * (N - M) / (M + 1);
}

double NegativeHyperGeometricRand::Variance() const
{
    double Mp1 = M + 1;
    double var = 1 - m / Mp1;
    var *= (N + 1) * (N - M);
    var /= Mp1 * (Mp1 + 1);
    return m * var;
}

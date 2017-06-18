#include "TrinomialRand.h"

TrinomialRand::TrinomialRand(int number, double probability1, double probability2)
{
    SetParameters(number, probability1, probability2);
}

std::string TrinomialRand::Name() const
{
    return "Trinomial(" + toStringWithPrecision(GetNumber()) + ", "
                        + toStringWithPrecision(GetFirstProbability()) + ", "
                        + toStringWithPrecision(GetSecondProbability()) + ")";
}

void TrinomialRand::SetParameters(int number, double probability1, double probability2)
{
    X.SetParameters(number, probability1);
    Y.SetParameters(number, probability2);

    n = X.GetNumber();
    log1mProb = std::log1p(-X.GetProbability() - Y.GetProbability());
    p2_1mp1 = Y.GetProbability() / (1.0 - X.GetProbability());
}

double TrinomialRand::P(const IntPair &point) const
{
    int x = point.first, y = point.second;
    return (x < 0 || x > n || y < 0 || y > n) ? 0.0 : std::exp(logP(point));
}

double TrinomialRand::logP(const IntPair &point) const
{
    int x = point.first, y = point.second;
    if (x < 0 || x > n || y < 0 || y > n)
        return -INFINITY;
    double res = X.GetLogFactorialN();
    res -= RandMath::lfact(x);
    res -= RandMath::lfact(y);
    res -= RandMath::lfact(n - x - y);
    res += x * X.GetLogProbability();
    res += y * Y.GetLogProbability();
    res += (n - x - y) * log1mProb;
    return res;
}

double TrinomialRand::F(const IntPair &point) const
{
    int x = point.first, y = point.second;
    if (x < 0 || y < 0)
        return 0.0;
    if (x >= n)
        return Y.F(y);
    if (y >= n)
        return X.F(x);
    double sum = 0.0;
    for (int i = 0; i <= x; ++i) {
        for (int j = 0; j <= y; ++j)
            sum += P(DoublePair(i, j));
    }
    return sum;
}

IntPair TrinomialRand::Variate() const
{
    int x = X.Variate();
    int y = BinomialRand::Variate(n - x, p2_1mp1);
    return IntPair(x, y);
}

double TrinomialRand::Correlation() const
{
    double p1 = X.GetProbability(), p2 = Y.GetProbability();
    double corr = p1 * p2;
    corr /= 1.0 - p1;
    corr /= 1.0 - p2;
    return -std::sqrt(corr);
}


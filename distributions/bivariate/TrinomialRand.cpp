#include "TrinomialRand.h"

TrinomialRand::TrinomialRand(int number, double probability1, double probability2)
{
    SetParameters(number, probability1, probability2);
}

String TrinomialRand::Name() const
{
    return "Trinomial(" + toStringWithPrecision(GetNumber()) + ", "
                        + toStringWithPrecision(GetFirstProbability()) + ", "
                        + toStringWithPrecision(GetSecondProbability()) + ")";
}

void TrinomialRand::SetParameters(int number, double probability1, double probability2)
{
    if (probability1 < 0.0 || probability2 < 0.0 || probability1 + probability2 > 1.0)
        throw std::invalid_argument("Probabilities of Trinomial distribution should be positive and their sum should be be not greater than 1");
    if (number <= 0)
        throw std::invalid_argument("Number of Trinomial distribution should be positive");

    double p1 = probability1, p2 = probability2;
    n = number;
    X.SetParameters(n, p1);
    Y.SetParameters(n, p2);

    log1mProb = std::log1p(-p1 - p2);
    p1_1mp2 = p1 / (1.0 - p2);
    p2_1mp1 = p2 / (1.0 - p1);
}

double TrinomialRand::P(const IntPair &point) const
{
    int x = point.first, y = point.second;
    return (x < 0 || x > n || y < 0 || y > n) ? 0.0 : std::exp(logP(point));
}

double TrinomialRand::logP(const IntPair &point) const
{
    int x = point.first, y = point.second;
    if (x < 0 || y < 0 || x + y > n)
        return -INFINITY;
    int nmxmy = n - x - y;
    double res = X.GetLogFactorialN();
    res -= RandMath::lfact(x);
    res -= RandMath::lfact(y);
    res -= RandMath::lfact(nmxmy);
    res += x * X.GetLogProbability();
    res += y * Y.GetLogProbability();
    res += nmxmy * log1mProb;
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
    double p1 = X.GetProbability(), p2 = Y.GetProbability();
    if (x > y) {
        double ratio = (1 - p2 - p1) / (1 - p2);
        for (int i = 0; i <= y; ++i) {
            int nmimx = n - i - x;
            double prob1 = (nmimx > 0) ? RandMath::ibeta(ratio, nmimx, x + 1) : 1.0;
            double prob2 = Y.P(i);
            sum += prob1 * prob2;
        }
        return sum;
    }
    double ratio = (1 - p1 - p2) / (1 - p1);
    for (int i = 0; i <= x; ++i) {
        int nmimy = n - i - y;
        double prob1 = (nmimy > 0) ? RandMath::ibeta(ratio, nmimy, y + 1) : 1.0;
        double prob2 = X.P(i);
        sum += prob1 * prob2;
    }
    return sum;
}

IntPair TrinomialRand::Variate() const
{
    if (X.GetProbability() > Y.GetProbability()) {
        int x = X.Variate();
        int y = BinomialRand::Variate(n - x, p2_1mp1);
        return IntPair(x, y);
    }
    int y = Y.Variate();
    int x = BinomialRand::Variate(n - y, p1_1mp2);
    return IntPair(x, y);
}

double TrinomialRand::Correlation() const
{
    double logp1 = X.GetLogProbability(), logp2 = Y.GetLogProbability();
    double log1mp1 = X.GetLog1mProbability(), log1mp2 = Y.GetLog1mProbability();
    return -std::exp(logp1 + logp2 - log1mp1 - log1mp2);
}


#include "TrinomialRand.h"

#include <algorithm>

template< typename IntType >
TrinomialRand<IntType>::TrinomialRand(int number, double probability1, double probability2)
{
    SetParameters(number, probability1, probability2);
}

template< typename IntType >
String TrinomialRand<IntType>::Name() const
{
    return "Trinomial(" + this->toStringWithPrecision(GetNumber()) + ", "
                        + this->toStringWithPrecision(GetFirstProbability()) + ", "
                        + this->toStringWithPrecision(GetSecondProbability()) + ")";
}

template< typename IntType >
void TrinomialRand<IntType>::SetParameters(int number, double probability1, double probability2)
{
    if (!(probability1 >= 0.0 && probability2 >= 0.0 && probability1 + probability2 <= 1.0))
        throw std::invalid_argument("Probabilities of Trinomial distribution should be positive and their sum should be not greater than 1, but they're equal to "
                                    + std::to_string(probability1) + " and " + std::to_string(probability2));
    if (!(number > 0))
        throw std::invalid_argument("Number of Trinomial distribution should be positive, but it's equal to "
                                    + std::to_string(number));

    double p1 = probability1, p2 = probability2;
    n = number;
    this->X.SetParameters(n, p1);
    this->Y.SetParameters(n, p2);

    log1mProb = std::log1pl(-p1 - p2);
    p1_1mp2 = (p2 == 1.0) ? 0.0 : p1 / (1.0 - p2);
    p2_1mp1 = (p1 == 1.0) ? 0.0 : p2 / (1.0 - p1);
}

template< typename IntType >
double TrinomialRand<IntType>::P(const Pair<IntType> &point) const
{
    IntType x = point.first, y = point.second;
    return (x < 0 || x > n || y < 0 || y > n) ? 0.0 : std::exp(logP(point));
}

template< typename IntType >
double TrinomialRand<IntType>::logP(const Pair<IntType> &point) const
{
    IntType x = point.first, y = point.second;
    if (x < 0 || y < 0 || x > n || y > n - x)
        return -INFINITY;
    IntType nmxmy = n - x - y;
    double res = this->X.GetLogFactorialN();
    res -= RandMath::lfact(x);
    res -= RandMath::lfact(y);
    res -= RandMath::lfact(nmxmy);
    if (x > 0)
        res += x * this->X.GetLogProbability();
    if (y > 0)
        res += y * this->Y.GetLogProbability();
    if (nmxmy > 0)
        res += nmxmy * log1mProb;
    return res;
}

template< typename IntType >
double TrinomialRand<IntType>::F(const Pair<IntType> &point) const
{
    int x = point.first, y = point.second;
    if (x < 0 || y < 0)
        return 0.0;
    double p1 = this->X.GetProbability(), p2 = this->Y.GetProbability();
    if (p1 == 1.0)
        return (x >= n) ? 1.0 : 0.0;
    if (p2 == 1.0)
        return (y >= n) ? 1.0 : 0.0;
    if (x >= n)
        return this->Y.F(y);
    if (y >= n)
        return this->X.F(x);
    double sum = 0.0;
    if (x > y) {
        double ratio = (1 - p2 - p1) / (1 - p2);
        for (int i = 0; i <= y; ++i) {
            int nmimx = n - i - x;
            double prob1 = (nmimx > 0) ? RandMath::ibeta(ratio, nmimx, x + 1) : 1.0;
            double prob2 = this->Y.P(i);
            sum += prob1 * prob2;
        }
        return sum;
    }
    double ratio = (1 - p1 - p2) / (1 - p1);
    for (int i = 0; i <= x; ++i) {
        IntType nmimy = n - i - y;
        double prob1 = (nmimy > 0) ? RandMath::ibeta(ratio, nmimy, y + 1) : 1.0;
        double prob2 = this->X.P(i);
        sum += prob1 * prob2;
    }
    return sum;
}

template< typename IntType >
Pair<IntType> TrinomialRand<IntType>::Variate() const
{
    if (this->X.GetProbability() == 1.0)
        return IntPair(n, 0);
    if (this->Y.GetProbability() == 1.0)
        return IntPair(0, n);
    if (this->X.GetProbability() > this->Y.GetProbability()) {
        IntType x = this->X.Variate();
        IntType y = BinomialDistribution<IntType>::Variate(n - x, p2_1mp1, this->localRandGenerator);
        return IntPair(x, y);
    }
    IntType y = this->Y.Variate();
    IntType x = BinomialDistribution<IntType>::Variate(n - y, p1_1mp2, this->localRandGenerator);
    return IntPair(x, y);
}

template< typename IntType >
long double TrinomialRand<IntType>::Correlation() const
{
    long double p1 = this->X.GetProbability(), p2 = this->Y.GetProbability();
    if (!(p1 > 0.0L && p1 < 1.0L && p2 > 0.0L && p2 < 1.0L))
        return NAN;
    return -std::sqrt(p1 * p2 / ((1.0L - p1) * (1.0L - p2)));
}

template< typename IntType >
Pair<IntType> TrinomialRand<IntType>::Mode() const
{
    IntType modeX = static_cast<IntType>(std::floor(n * this->X.GetProbability()));
    IntType modeY = static_cast<IntType>(std::floor(n * this->Y.GetProbability()));
    modeY = std::min(modeY, static_cast<IntType>(n - modeX));
    Pair<IntType> mode(modeX, modeY);
    double modeValue = logP(mode);

    while (true) {
        Pair<IntType> best = mode;
        double bestValue = modeValue;
        auto consider = [this, &best, &bestValue](IntType x, IntType y) {
            double value = logP(std::make_pair(x, y));
            if (value > bestValue) {
                best = std::make_pair(x, y);
                bestValue = value;
            }
        };

        if (mode.first + mode.second < n) {
            consider(mode.first + 1, mode.second);
            consider(mode.first, mode.second + 1);
        }
        if (mode.first > 0) {
            consider(mode.first - 1, mode.second);
            consider(mode.first - 1, mode.second + 1);
        }
        if (mode.second > 0) {
            consider(mode.first, mode.second - 1);
            consider(mode.first + 1, mode.second - 1);
        }
        if (best == mode)
            return mode;
        mode = best;
        modeValue = bestValue;
    }
}

template class TrinomialRand<int>;
template class TrinomialRand<long int>;
template class TrinomialRand<long long int>;


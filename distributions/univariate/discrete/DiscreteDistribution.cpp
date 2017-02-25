#include "DiscreteDistribution.h"

void DiscreteDistribution::ProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

int DiscreteDistribution::Mode() const
{
    /// Works only for unimodal and monotone from start point to mode distributions
    int x = Median();
    double prob = P(x), newProb = P(x + 1);
    if (prob < newProb) {
        do {
            ++x;
            prob = newProb;
            newProb = P(x + 1);
        } while (prob < newProb);
    }
    else {
        newProb = P(x - 1);
        while (prob < newProb) {
            --x;
            prob = newProb;
            newProb = P(x - 1);
        }
    }
    return x;
}

double DiscreteDistribution::quantileImpl(double p) const
{
    // TODO: make sample quantile as starting point
    double mean = Mean();
    int down = static_cast<int>(std::floor(mean)), up = down + 1;
    double fu = F(up), fd = F(down);
    int iter = 0;
    static constexpr int maxIter = 1e4;
    /// go up
    while (fu < p)
    {
        fd = fu;
        fu = F(++up);
        if (++iter > maxIter)
            return NAN;
    }
    down = up - 1;

    iter = 0;
    /// go down
    while (fd > p)
    {
        fd = F(--down);
        if (++iter > maxIter)
            return NAN;
    }
    up = down + 1;

    /// if lower quantile is not equal probability, we return upper quantile
    return (fd < p) ? up : down;
}

double DiscreteDistribution::quantileImpl1m(double p) const
{
    double mean = Mean();
    int down = static_cast<int>(std::floor(mean)), up = down + 1;
    double su = S(up), sd = S(down);
    int iter = 0;
    static constexpr int maxIter = 1e4;
    /// go up
    while (su > p)
    {
        sd = su;
        su = S(++up);
        if (++iter > maxIter)
            return NAN;
    }
    down = up - 1;

    iter = 0;
    /// go down
    while (sd < p)
    {
        sd = S(--down);
        if (++iter > maxIter)
            return NAN;
    }
    up = down + 1;

    /// if lower quantile is not equal probability, we return upper quantile
    return (sd > p) ? up : down;
}

double DiscreteDistribution::ExpectedValue(const std::function<double (double)> &funPtr, int minPoint, int maxPoint) const
{
    SUPPORT_TYPE suppType = SupportType();
    int k = minPoint, upperBoundary = maxPoint;
    if (suppType == FINITE_T || suppType == RIGHTSEMIFINITE_T) {
        k = std::max(k, static_cast<int>(std::floor(MinValue())));
    }
    else if (k == INT_MIN) { /// neither function nor distribution are bounded from the left
        k = Quantile(1e-6);
    }
    if (suppType == FINITE_T || suppType == LEFTSEMIFINITE_T) {
        upperBoundary = std::min(upperBoundary, static_cast<int>(std::ceil(MaxValue())));
    }
    else if (upperBoundary == INT_MAX) { /// neither function nor distribution are bounded from the right
        upperBoundary = Quantile1m(1e-6);
    }

    double sum = 0;
    do {
        double addon = funPtr(k);
        if (addon != 0.0) {
            addon *= P(k);
            sum += addon;
        }
        ++k;
    } while (k <= upperBoundary);
    return sum;
}

double DiscreteDistribution::Hazard(double x) const
{
    if (x < MinValue())
        return 0.0; /// 0/1
    if (x > MaxValue())
        return NAN; /// 0/0
    return P(x) / S(x);
}

double DiscreteDistribution::Likelihood(const std::vector<int> &sample) const
{
    double res = 1.0;
    for (const int & var : sample )
        res *= P(var);
    return res;
}

double DiscreteDistribution::LogLikelihood(const std::vector<int> &sample) const
{
    double res = 0.0;
    for (const int & var : sample )
        res += logP(var);
    return res;
}

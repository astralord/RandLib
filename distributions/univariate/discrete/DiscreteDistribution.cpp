#include "DiscreteDistribution.h"
#include "../continuous/GammaRand.h"

void DiscreteDistribution::ProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

void DiscreteDistribution::LogProbabilityMassFunction(const std::vector<int> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = logP(x[i]);
}

int DiscreteDistribution::Mode() const
{
    /// Works only for unimodal and monotone from starting point to the mode distributions
    int x = Median();
    double logProb = logP(x), newLogProb = logP(x + 1);
    if (logProb < newLogProb) {
        do {
            ++x;
            logProb = newLogProb;
            newLogProb = logP(x + 1);
        } while (logProb < newLogProb);
    }
    else {
        newLogProb = logP(x - 1);
        while (logProb < newLogProb) {
            --x;
            logProb = newLogProb;
            newLogProb = logP(x - 1);
        }
    }
    return x;
}

int DiscreteDistribution::quantileImpl(double p) const
{
    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<int> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    std::nth_element(sample.begin(), sample.begin() + index, sample.end());
    int guess = sample[index];
    int down = static_cast<int>(std::floor(guess)), up = down + 1;
    double fu = F(up), fd = F(down);
    /// go up
    while (fu < p) {
        fd = fu;
        fu = F(++up);
    }
    down = up - 1;
    /// go down
    while (fd > p) {
        fd = F(--down);
    }
    up = down + 1;
    /// if lower quantile is not equal probability, we return upper quantile
    return (fd < p) ? up : down;
}

int DiscreteDistribution::quantileImpl1m(double p) const
{
    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<int> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    std::nth_element(sample.begin(), sample.begin() + index, sample.end(), std::greater<>());
    int guess = sample[index];
    int down = static_cast<int>(std::floor(guess)), up = down + 1;
    double su = S(up), sd = S(down);
    /// go up
    while (su > p) {
        sd = su;
        su = S(++up);
    }
    down = up - 1;
    /// go down
    while (sd < p) {
        sd = S(--down);
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

double DiscreteDistribution::LikelihoodFunction(const std::vector<int> &sample) const
{
    double res = 1.0;
    for (const int & var : sample )
        res *= P(var);
    return res;
}

double DiscreteDistribution::LogLikelihoodFunction(const std::vector<int> &sample) const
{
    double res = 0.0;
    for (const int & var : sample )
        res += logP(var);
    return res;
}

bool DiscreteDistribution::PearsonChiSquaredTest(const std::vector<int> &orderStatistic, double alpha, int lowerBoundary, int upperBoundary, size_t numberOfEstimatedParameters) const
{
    size_t n = orderStatistic.size(), i = 0, k = 0;
    double nInv = 1.0 / n, sum = 0.0;

    /// Sanity checks
    if (lowerBoundary >= upperBoundary)
        throw std::invalid_argument("Lower boundary should be less than upper one");
    for (size_t j = 1; j != n; ++j) {
        if (orderStatistic[i] < orderStatistic[j - 1])
            throw std::invalid_argument("Sample should be sorted in ascending order");
    }
    if (orderStatistic[0] < this->MinValue())
        throw std::invalid_argument("Some elements in the sample are too small to belong to this distribution, they should be bigger than " + toStringWithPrecision(this->MinValue()));
    if (orderStatistic[n - 1] > this->MaxValue())
        throw std::invalid_argument("Some elements in the sample are too large to belong to this distribution, they should be less than " + toStringWithPrecision(this->MaxValue()));

    /// Lower interval
    int x = orderStatistic[0];
    if (lowerBoundary > this->MinValue()) {
        auto upIt = std::upper_bound(orderStatistic.begin(), orderStatistic.end(), lowerBoundary);
        i += upIt - orderStatistic.begin();
        x = orderStatistic[i];
        double prob = nInv * i, expectedProb = this->F(lowerBoundary);
        double addon = prob - expectedProb;
        addon *= addon;
        addon /= expectedProb;
        sum += addon;
        ++k;
    }
    /// Middle intervals
    while (i < n && x < upperBoundary) {
        size_t count = 1;
        x = orderStatistic[i];
        while (i + count < n && x == orderStatistic[i + count])
            ++count;
        double prob = nInv * count, expectedProb = this->P(x);
        double addon = prob - expectedProb;
        addon *= addon;
        addon /= expectedProb;
        sum += addon;
        i += count;
        ++k;
    }
    /// Upper interval
    if (upperBoundary < this->MaxValue()) {
        double prob = nInv * (n - i), expectedProb = this->S(upperBoundary);
        double addon = prob - expectedProb;
        addon *= addon;
        addon /= expectedProb;
        sum += addon;
        ++k;
    }

    if (k <= numberOfEstimatedParameters + 1) {
        throw std::invalid_argument("Sample is too small, number of groups (" + toStringWithPrecision(k)
                                    + ") should be bigger than number of estimated parameters plus one (" + toStringWithPrecision(numberOfEstimatedParameters + 1) + ")");
    }
    double statistic = n * sum;
    ChiSquaredRand X(k - 1);
    double q = X.Quantile1m(alpha);
    return (statistic <= q);
}

bool DiscreteDistribution::PearsonChiSquaredTest(const std::vector<int> &orderStatistic, double alpha, size_t numberOfEstimatedParameters) const
{
    return PearsonChiSquaredTest(orderStatistic, alpha, this->MinValue(), this->MaxValue(), numberOfEstimatedParameters);
}

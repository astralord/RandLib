#include "DiscreteDistribution.h"
#include "../continuous/GammaRand.h"

template < typename IntType >
void DiscreteDistribution<IntType>::ProbabilityMassFunction(const std::vector<IntType> &x, std::vector<double> &y) const
{
    for (size_t i = 0; i != x.size(); ++i)
        y[i] = this->P(x[i]);
}

template < typename IntType >
void DiscreteDistribution<IntType>::LogProbabilityMassFunction(const std::vector<IntType> &x, std::vector<double> &y) const
{
    for (size_t i = 0; i != x.size(); ++i)
        y[i] = this->logP(x[i]);
}

template < typename IntType >
IntType DiscreteDistribution<IntType>::Mode() const
{
    /// Works only for unimodal distributions
    IntType x = this->Median();
    double logProb = this->logP(x), newLogProb = this->logP(x + 1);
    if (logProb < newLogProb) {
        do {
            ++x;
            logProb = newLogProb;
            newLogProb = this->logP(x + 1);
        } while (logProb < newLogProb);
    }
    else {
        newLogProb = this->logP(x - 1);
        while (logProb < newLogProb) {
            --x;
            logProb = newLogProb;
            newLogProb = this->logP(x - 1);
        }
    }
    return x;
}

template < typename IntType >
IntType DiscreteDistribution<IntType>::quantileImpl(double p, IntType initValue) const
{
    IntType down = initValue, up = down + 1;
    double fu = this->F(up), fd = this->F(down);
    /// go up
    while (fu < p) {
        fd = fu;
        fu = this->F(++up);
    }
    down = up - 1;
    /// go down
    while (fd > p) {
        fd = this->F(--down);
    }
    up = down + 1;
    /// if lower quantile is not equal probability, we return upper quantile
    return (fd < p) ? up : down;
}


template < typename IntType >
IntType DiscreteDistribution<IntType>::quantileImpl(double p) const
{
    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<IntType> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    if (index == 0)
        return this->quantileImpl(p, *std::min_element(sample.begin(), sample.end()));
    std::nth_element(sample.begin(), sample.begin() + index, sample.end());
    return this->quantileImpl(p, sample[index]);
}

template < typename IntType >
IntType DiscreteDistribution<IntType>::quantileImpl1m(double p, IntType initValue) const
{
    IntType down = initValue, up = down + 1;
    double su = this->S(up), sd = this->S(down);
    /// go up
    while (su > p) {
        sd = su;
        su = this->S(++up);
    }
    down = up - 1;
    /// go down
    while (sd < p) {
        sd = this->S(--down);
    }
    up = down + 1;

    /// if lower quantile is not equal probability, we return upper quantile
    return (sd > p) ? up : down;
}

template < typename IntType >
IntType DiscreteDistribution<IntType>::quantileImpl1m(double p) const
{
    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<IntType> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    if (index == 0)
        return this->quantileImpl1m(p, *std::max_element(sample.begin(), sample.end()));
    std::nth_element(sample.begin(), sample.begin() + index, sample.end(), std::greater<>());
    return this->quantileImpl1m(p, sample[index]);
}

template < typename IntType >
long double DiscreteDistribution<IntType>::ExpectedValue(const std::function<double (IntType)> &funPtr, IntType minPoint, IntType maxPoint) const
{
    SUPPORT_TYPE suppType = this->SupportType();
    IntType k = minPoint, upperBoundary = maxPoint;
    if (suppType == FINITE_T || suppType == RIGHTSEMIFINITE_T) {
        k = std::max(k, this->MinValue());
    }
    if (suppType == FINITE_T || suppType == LEFTSEMIFINITE_T) {
        upperBoundary = std::min(upperBoundary, this->MaxValue());
    }

    double sum = 0;
    do {
        double addon = funPtr(k);
        if (addon != 0.0) {
            double prob = this->P(k);
            if (prob < MIN_POSITIVE)
                return sum;
            addon *= this->P(k);
            sum += addon;
        }
        ++k;
    } while (k <= upperBoundary);
    return sum;
}

template < typename IntType >
double DiscreteDistribution<IntType>::Hazard(const IntType &x) const
{
    if (x < this->MinValue())
        return 0.0; /// 0/1
    if (x > this->MaxValue())
        return NAN; /// 0/0
    return this->P(x) / this->S(x);
}

template < typename IntType >
double DiscreteDistribution<IntType>::LikelihoodFunction(const std::vector<IntType> &sample) const
{
    long double res = 1.0;
    for (const IntType & var : sample )
        res *= this->P(var);
    return res;
}

template < typename IntType >
double DiscreteDistribution<IntType>::LogLikelihoodFunction(const std::vector<IntType> &sample) const
{
    long double res = 0.0;
    for (const IntType & var : sample )
        res += this->logP(var);
    return res;
}

template < typename IntType >
bool DiscreteDistribution<IntType>::PearsonChiSquaredTest(const std::vector<IntType> &orderStatistic, double alpha, int lowerBoundary, int upperBoundary, size_t numberOfEstimatedParameters) const
{
    size_t n = orderStatistic.size(), i = 0, k = 0;
    double nInv = 1.0 / n, sum = 0.0;

    /// Sanity checks
    if (lowerBoundary >= upperBoundary)
        throw std::invalid_argument("Lower boundary should be smaller than upper one");
    for (size_t j = 1; j != n; ++j) {
        if (orderStatistic[i] < orderStatistic[j - 1])
            throw std::invalid_argument("Order statistic should be sorted in ascending order");
    }
    if (orderStatistic[0] < this->MinValue())
        throw std::invalid_argument("Some elements in the sample are too small to belong to this distribution, they should be larger than " + this->toStringWithPrecision(this->MinValue()));
    if (orderStatistic[n - 1] > this->MaxValue())
        throw std::invalid_argument("Some elements in the sample are too large to belong to this distribution, they should be smaller than " + this->toStringWithPrecision(this->MaxValue()));

    /// Lower interval
    IntType x = orderStatistic[0];
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
        throw std::invalid_argument("Sample is too small, number of groups (" + this->toStringWithPrecision(k)
                                    + ") should be bigger than number of estimated parameters plus 1 ("
                                    + this->toStringWithPrecision(numberOfEstimatedParameters + 1) + ")");
    }
    double statistic = n * sum;
    ChiSquaredRand X(k - 1);
    double q = X.Quantile1m(alpha);
    return statistic <= q;
}

template < typename IntType >
bool DiscreteDistribution<IntType>::PearsonChiSquaredTest(const std::vector<IntType> &orderStatistic, double alpha, size_t numberOfEstimatedParameters) const
{
    return PearsonChiSquaredTest(orderStatistic, alpha, this->MinValue(), this->MaxValue(), numberOfEstimatedParameters);
}

template class DiscreteDistribution<int>;
template class DiscreteDistribution<long int>;
template class DiscreteDistribution<long long int>;

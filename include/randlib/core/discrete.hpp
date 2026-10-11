#pragma once

#include <span>

#include <limits>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/univariate.hpp>

namespace randlib
{

/**
 *@brief The DiscreteDistribution class <BR>
 * Abstract class for all discrete distributions
 */
template <typename IntType>
class DiscreteDistribution : virtual public UnivariateDistribution<IntType>
{
    static_assert(std::is_integral_v<IntType> && std::is_signed_v<IntType>, "Discrete distribution supports only signed integral types");

protected:
    DiscreteDistribution() {}
public:
    virtual ~DiscreteDistribution() {}
protected:

public:
    /**
     * @fn P
     * @param k
     * @return probability to get k
     */
    [[nodiscard]] virtual double P(const IntType & k) const { return std::exp(this->logP(k)); }

    /**
     * @fn logP
     * @param x
     * @return logarithm of probability to get x
     */
    [[nodiscard]] virtual double logP(const IntType & x) const = 0;

    /**
     * @fn ProbabilityMassFunction
     * fill vector y with P(x)
     * @param x
     * @param y
     */
    void ProbabilityMassFunction(std::span<const IntType> x, std::span<double> y) const;

    /**
     * @fn LogProbabilityMassFunction
     * fill vector y with logP(x)
     * @param x
     * @param y
     */
    void LogProbabilityMassFunction(std::span<const IntType> x, std::span<double> y) const;

    [[nodiscard]] IntType Mode() const override;

protected:
    IntType quantileImpl(double p, IntType initValue) const override;
    IntType quantileImpl(double p) const override;
    IntType quantileImpl1m(double p, IntType initValue) const override;
    IntType quantileImpl1m(double p) const override;
    long double ExpectedValue(const std::function<double (IntType)> &funPtr, IntType minPoint, IntType maxPoint) const override;
    std::complex<double> CFImpl(double t) const override;

private:
    /**
     * @fn searchQuantile
     * Deterministic exponential and binary search over the support
     * @param p probability in (0, 1)
     * @param upperTail search for the smallest k with S(k) ≤ p instead of F(k) ≥ p
     * @return the smallest k with F(k) ≥ p (or S(k) ≤ p)
     */
    IntType searchQuantile(double p, bool upperTail) const;

public:
    /**
     * @fn Hazard
     * @param x
     * @return hazard function
     */
    [[nodiscard]] double Hazard(const IntType &x) const override;

    /**
     * @fn LikelihoodFunction
     * @param sample
     * @return likelihood function of the distribution for given sample
     */
    [[nodiscard]] double LikelihoodFunction(std::span<const IntType> sample) const override;

    /**
     * @fn LogLikelihoodFunction
     * @param sample
     * @return log-likelihood function of the distribution for given sample
     */
    [[nodiscard]] double LogLikelihoodFunction(std::span<const IntType> sample) const override;

    /**
     * @fn PearsonChiSquaredTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha significance level of the test
     * @param lowerBoundary setting of the left most interval (-infinity, lowerBoundary]
     * @param upperBoundary setting of the right most interval [upperBoundary, infinity)
     * @param numberOfEstimatedParameters zero by default
     * @return true if sample is from this distribution according to Pearson's chi-squared test, false otherwise
     *
     * Chi-squared test is implemented only for discrete distribution, as it is much more complicated for continuous one.
     * The problem is ambiguity of grouping sample into intervals. In the case when parameters are estimated by
     * maximum-likelihood estimator, using original observations, statistics might not follow asymptotic chi-square
     * distribution and that leads to serious underestimate of the error of the first kind.
     * For more details look: "The use of MLE in chi-square tests for goodness of fit" by Herman Chernoff and E.L. Lehmann
     */
    bool PearsonChiSquaredTest(std::span<const IntType> orderStatistic, double alpha, int lowerBoundary, int upperBoundary, size_t numberOfEstimatedParameters = 0) const;

    /**
     * @fn PearsonChiSquaredTest
     * @param orderStatistic sample sorted in ascending order
     * @param alpha significance level of the test
     * @param numberOfEstimatedParameters zero by default
     * @return true if sample is from this distribution according to Pearson's chi-squared test, false otherwise
     * In this function user won't set upper and lower intervals for tails.
     * However it might be useful to group rare events for chi-squared test to achieve better results
     */
    bool PearsonChiSquaredTest(std::span<const IntType> orderStatistic, double alpha, size_t numberOfEstimatedParameters = 0) const;
};

} // namespace randlib

#include <randlib/core/exponential_family.hpp>
#include <randlib/continuous/gamma.hpp>

namespace randlib
{

template < typename IntType >
void DiscreteDistribution<IntType>::ProbabilityMassFunction(std::span<const IntType> x, std::span<double> y) const
{
    for (size_t i = 0; i != x.size(); ++i)
        y[i] = this->P(x[i]);
}

template < typename IntType >
void DiscreteDistribution<IntType>::LogProbabilityMassFunction(std::span<const IntType> x, std::span<double> y) const
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
  return searchQuantile(p, false);
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
  return searchQuantile(p, true);
}

template < typename IntType >
IntType DiscreteDistribution<IntType>::searchQuantile(double p, bool upperTail) const
{
    /// reached(k) is false below the quantile and true from it on
    auto reached = [this, p, upperTail] (long long k) {
      return upperTail ? this->S(static_cast<IntType>(k)) <= p : this->F(static_cast<IntType>(k)) >= p;
    };
    const long long lo = this->MinValue(), hi = this->MaxValue();
    const bool loBounded = lo > static_cast<long long>(std::numeric_limits<IntType>::lowest());

    /// Find a < b with !reached(a) and reached(b) by exponential search from the bounded end or 0
    long long a, b;
    if (loBounded || !reached(0)) {
        a = loBounded ? lo - 1 : 0;
        long long step = 1;
        b = a + step;
        while (b < hi && !reached(b)) {
            a = b;
            step *= 2;
            b = (step > hi - a) ? hi : a + step;
        }
        b = std::min(b, hi);
    }
    else {
        b = 0;
        long long step = 1;
        a = b - step;
        while (a > lo && reached(a)) {
            b = a;
            step *= 2;
            a = (step > a - lo) ? lo : b - step;
        }
    }
    while (b - a > 1) {
        long long middle = a + (b - a) / 2;
        if (reached(middle))
            b = middle;
        else
            a = middle;
    }
  return static_cast<IntType>(b);
}

template < typename IntType >
long double DiscreteDistribution<IntType>::ExpectedValue(const std::function<double (IntType)> &funPtr, IntType minPoint, IntType maxPoint) const
{
    Support suppType = this->SupportType();
    IntType k = minPoint, upperBoundary = maxPoint;
    if (suppType == Support::Finite || suppType == Support::RightSemiFinite) {
        k = std::max(k, this->MinValue());
    }
    if (suppType == Support::Finite || suppType == Support::LeftSemiFinite) {
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
std::complex<double> DiscreteDistribution<IntType>::CFImpl(double t) const
{
    /// Σ P(k) e^{itk}, summed outwards from the median. Beyond the last term K of a walk with
    /// ratio r = e^{±it}, the remainder Σ_{m≥1} P(K + m) r^{K+m} is replaced by the first terms of
    /// Euler's transformation r^{K+1} / (1 - r) Σ_j (r / (1 - r))^j Δ^j P(K + 1), with forward
    /// differences Δ along the walk. The walk stops once the next term of that series is negligible
    /// or the remaining mass is; summing until P(k) itself vanishes would take ~10^10 terms for Zeta(2)
    static constexpr long long MAX_TERMS = 1000000;
    static constexpr double EPSILON = 1e-17;
    const long long lo = this->MinValue(), hi = this->MaxValue(), start = this->Median();
    std::complex<long double> sum = 0;
    auto walk = [&] (long long from, long long step, long long end, bool upper) {
        const std::complex<double> ratio = std::polar(1.0, step * t);
        const std::complex<double> w = ratio / (1.0 - ratio);
        const double damping = std::abs(1.0 - ratio);
        std::complex<double> zk = std::polar(1.0, t * static_cast<double>(from));
        double last[4] = {0, 0, 0, 0};
        long long k = from;
        for (long long terms = 0; (step > 0) ? k <= end : k >= end; k += step, ++terms) {
            double p = this->P(static_cast<IntType>(k));
            sum += std::complex<long double>(p * zk);
            zk *= ratio;
            if (terms % 64 == 63)
                zk = std::polar(1.0, t * static_cast<double>(k + step));
            last[3] = last[2]; last[2] = last[1]; last[1] = last[0]; last[0] = p;
            if (p == 0.0 && terms > 16)
              return;
            if (terms < 16 || damping == 0.0)
                continue;
            double thirdDifference = std::fabs(last[0] - 3 * last[1] + 3 * last[2] - last[3]);
            double nextTerm = std::pow(std::abs(w), 3) * thirdDifference / damping;
            bool unbounded = (step > 0) ? end - k > MAX_TERMS : k - end > MAX_TERMS;
            bool decreasing = last[0] < last[1] && last[1] < last[2] && last[2] < last[3];
            bool converged = unbounded && decreasing && nextTerm < EPSILON * std::abs(sum);
            bool powerOfTwo = ((terms + 1) & terms) == 0;
            bool massGone = powerOfTwo
                    && (upper ? this->S(static_cast<IntType>(k)) : this->F(static_cast<IntType>(k - 1))) < EPSILON;
            if (converged || massGone || terms >= MAX_TERMS)
                break;
        }
        if ((step > 0) ? k > end : k < end)
          return;
        k += step;
        double next[3];
        for (int j = 0; j < 3; ++j)
            next[j] = ((step > 0) ? k + j * step <= end : k + j * step >= end) ? this->P(static_cast<IntType>(k + j * step)) : 0.0;
        std::complex<double> series = next[0] + w * (next[1] - next[0]) + w * w * (next[2] - 2 * next[1] + next[0]);
        sum += std::complex<long double>(zk * series / (1.0 - ratio));
    };
    walk(start, 1, hi, true);
    if (start > lo)
        walk(start - 1, -1, lo, false);
  return std::complex<double>(sum);
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
double DiscreteDistribution<IntType>::LikelihoodFunction(std::span<const IntType> sample) const
{
    long double res = 1.0;
    for (const IntType & var : sample )
        res *= this->P(var);
  return res;
}

template < typename IntType >
double DiscreteDistribution<IntType>::LogLikelihoodFunction(std::span<const IntType> sample) const
{
    long double res = 0.0;
    for (const IntType & var : sample )
        res += this->logP(var);
  return res;
}

template < typename IntType >
bool DiscreteDistribution<IntType>::PearsonChiSquaredTest(std::span<const IntType> orderStatistic, double alpha, int lowerBoundary, int upperBoundary, size_t numberOfEstimatedParameters) const
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
        throw std::invalid_argument("Some elements in the sample are too small to belong to this distribution, they should be greater than "
                                    + this->toStringWithPrecision(this->MinValue()));
    if (orderStatistic[n - 1] > this->MaxValue())
        throw std::invalid_argument("Some elements in the sample are too large to belong to this distribution, they should be smaller than "
                                    + this->toStringWithPrecision(this->MaxValue()));

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
bool DiscreteDistribution<IntType>::PearsonChiSquaredTest(std::span<const IntType> orderStatistic, double alpha, size_t numberOfEstimatedParameters) const
{
  return PearsonChiSquaredTest(orderStatistic, alpha, this->MinValue(), this->MaxValue(), numberOfEstimatedParameters);
}

} // namespace randlib

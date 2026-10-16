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
 * @brief The ContinuousDistribution class <BR>
 * Abstract class for all continuous distributions
 */
template <typename RealType>
class ContinuousDistribution : virtual public UnivariateDistribution<RealType>
{
    static_assert(std::is_floating_point_v<RealType>, "Continuous distribution supports only floating-point types");

protected:
    ContinuousDistribution() {}

public:
    virtual ~ContinuousDistribution() {}

protected:
public:
    /**
     * @param x
     * @return probability density function
     */
    [[nodiscard]] virtual double f(const RealType &x) const { return std::exp(this->logf(x)); }

    /**
     * @param x
     * @return logarithm of probability density function
     */
    [[nodiscard]] virtual double logf(const RealType &x) const = 0;

    /**
     * fill vector y with f(x)
     * @param x
     * @param y
     */
    void ProbabilityDensityFunction(std::span<const RealType> x, std::span<double> y) const;

    /**
     * fill vector y with logf(x)
     * @param x
     * @param y
     */
    void LogProbabilityDensityFunction(std::span<const RealType> x, std::span<double> y) const;

    [[nodiscard]] RealType Mode() const override;

protected:
    RealType quantileImpl(double p, RealType initValue) const override;
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p, RealType initValue) const override;
    RealType quantileImpl1m(double p) const override;
    long double ExpectedValue(const std::function<double(RealType)> &funPtr, RealType minPoint,
                              RealType maxPoint) const override;

    /**
     * Deterministic root search for the quantile: brackets the root, starting from the
     * guess when it is finite, then refines it with the Illinois variant of regula falsi
     * (on a logarithmic scale for small p), falling back to bisection whenever the bracket
     * does not halve. It stops when the bracket cannot shrink in floating point.
     * @param p probability in (0, 1)
     * @param upperTail search for S(x) = p instead of F(x) = p
     * @param guess initial guess, or NaN
     * @param relativeStep first step of the bracket expansion, relative to the guess
     * @return such x that F(x) = p (or S(x) = p)
     */
    RealType searchQuantile(double p, bool upperTail, RealType guess, double relativeStep = 1e-3) const;

public:
    [[nodiscard]] double Hazard(const RealType &x) const override;
    [[nodiscard]] double LikelihoodFunction(std::span<const RealType> sample) const override;
    [[nodiscard]] double LogLikelihoodFunction(std::span<const RealType> sample) const override;

    /**
     * @param orderStatistic sample sorted in ascending order
     * @param alpha level of test
     * @return true if sample is from this distribution according to asymptotic KS-test, false otherwise
     */
    bool KolmogorovSmirnovTest(std::span<const RealType> orderStatistic, double alpha) const;
};

} // namespace randlib

#include <randlib/core/exponential_family.hpp>
#include <randlib/continuous/kolmogorov_smirnov.hpp>

namespace randlib
{

template <typename RealType>
void ContinuousDistribution<RealType>::ProbabilityDensityFunction(std::span<const RealType> x,
                                                                  std::span<double> y) const
{
    for (size_t i = 0; i != x.size(); ++i)
        y[i] = this->f(x[i]);
}

template <typename RealType>
void ContinuousDistribution<RealType>::LogProbabilityDensityFunction(std::span<const RealType> x,
                                                                     std::span<double> y) const
{
    for (size_t i = 0; i != x.size(); ++i)
        y[i] = this->logf(x[i]);
}

template <typename RealType>
RealType ContinuousDistribution<RealType>::quantileImpl(double p, RealType initValue) const
{
    return searchQuantile(p, false, initValue);
}

template <typename RealType>
RealType ContinuousDistribution<RealType>::quantileImpl(double p) const
{
    return searchQuantile(p, false, std::numeric_limits<RealType>::quiet_NaN());
}

template <typename RealType>
RealType ContinuousDistribution<RealType>::quantileImpl1m(double p, RealType initValue) const
{
    return searchQuantile(p, true, initValue);
}

template <typename RealType>
RealType ContinuousDistribution<RealType>::quantileImpl1m(double p) const
{
    return searchQuantile(p, true, std::numeric_limits<RealType>::quiet_NaN());
}

template <typename RealType>
RealType ContinuousDistribution<RealType>::searchQuantile(double p, bool upperTail, RealType guess,
                                                          double relativeStep) const
{
    /// Solve in the tail with the smaller probability: 1 - p is exact for p ≥ 1/2, and
    /// F (or S) cannot resolve a probability close to 1 as well as the other tail resolves 1 - p
    if (p > 0.5) {
        p = 1.0 - p;
        upperTail = !upperTail;
    }
    /// g increases through zero at the quantile; for small p it is compared on a log scale
    const bool logScale = p < 1e-3;
    const double logP = std::log(p);
    auto g = [this, p, upperTail, logScale, logP](RealType x) {
        double v = upperTail ? this->S(x) : this->F(x);
        if (!logScale)
            return upperTail ? p - v : v - p;
        double logV = (v > 0) ? std::log(v) : -INFINITY;
        return upperTail ? logP - logV : logV - logP;
    };
    const RealType lo = this->MinValue(), hi = this->MaxValue();
    const bool loFinite = std::isfinite(lo), hiFinite = std::isfinite(hi);

    /// Bracket [a, b] with g(a) < 0 <= g(b), expanding geometrically from the guess,
    /// from a finite end of the support, or from 0
    const bool useGuess = std::isfinite(guess) && guess > lo && guess < hi;
    const RealType origin = useGuess ? guess : (loFinite ? lo : (hiFinite ? hi : 0));
    RealType a = lo, b = hi;
    double ga = -INFINITY, gb = INFINITY;
    if (useGuess || !(loFinite && hiFinite)) {
        double gOrigin = g(origin);
        if (gOrigin == 0)
            return origin;
        RealType step = useGuess ? std::max<RealType>(std::fabs(guess) * relativeStep, 1e-3) : 1;
        bool searchUp = (useGuess || !loFinite || hiFinite) ? gOrigin < 0 : true;
        if (!useGuess && loFinite && !hiFinite)
            searchUp = true;
        if (searchUp) {
            a = origin;
            ga = gOrigin;
            for (;;) {
                b = origin + step;
                if (hiFinite && b >= hi) {
                    b = hi;
                    gb = INFINITY;
                    break;
                }
                if (!std::isfinite(b))
                    return b;
                gb = g(b);
                if (gb >= 0)
                    break;
                a = b;
                ga = gb;
                step += step;
            }
        } else {
            b = origin;
            gb = gOrigin;
            for (;;) {
                a = origin - step;
                if (loFinite && a <= lo) {
                    a = lo;
                    ga = -INFINITY;
                    break;
                }
                if (!std::isfinite(a))
                    return a;
                ga = g(a);
                if (ga < 0)
                    break;
                b = a;
                gb = ga;
                step += step;
            }
        }
    }

    /// Bisection point: geometric in the distance to a finite support end while the
    /// bracket spans orders of magnitude, arithmetic otherwise
    const RealType anchor = loFinite ? lo : (hiFinite ? hi : 0);
    auto split = [anchor, loFinite, hiFinite](RealType left, RealType right) {
        const RealType direction = (loFinite || (!hiFinite && left >= 0)) ? 1 : -1;
        RealType tl = direction * (left - anchor), tr = direction * (right - anchor);
        RealType tNear = std::min(tl, tr), tFar = std::max(tl, tr), t;
        if (tNear > 0 && tFar > 4 * tNear)
            t = std::sqrt(tNear) * std::sqrt(tFar);
        else if (tNear == 0 && tFar > 1e-280 && (loFinite || hiFinite))
            t = tFar / 65536;
        else
            t = tNear + (tFar - tNear) / 2;
        return anchor + direction * t;
    };

    /// While the bracket spans orders of magnitude or g is infinite at an end, split it
    /// geometrically; then TOMS 748 refines it to full precision
    for (int iter = 0; iter < 2000; ++iter) {
        RealType tl = std::fabs(a - anchor), tr = std::fabs(b - anchor);
        bool spansMagnitudes =
            std::min(tl, tr) == 0 ? std::max(tl, tr) > 1e-280 : std::max(tl, tr) > 4 * std::min(tl, tr);
        if (std::isfinite(ga) && std::isfinite(gb) && !(spansMagnitudes && (loFinite || hiFinite)))
            break;
        RealType middle = split(a, b);
        if (!(middle > a && middle < b))
            return middle;
        double gm = g(middle);
        if (gm == 0)
            return middle;
        if (gm < 0) {
            a = middle;
            ga = gm;
        } else {
            b = middle;
            gb = gm;
        }
    }
    auto [left, right] = math::toms748<RealType>([&g](RealType x) { return static_cast<RealType>(g(x)); }, a, b,
                                                 static_cast<RealType>(ga), static_cast<RealType>(gb));
    return left + (right - left) / 2;
}

template <typename RealType>
RealType ContinuousDistribution<RealType>::Mode() const
{
    RealType guess = this->Mean(); /// good starting point
    if (!std::isfinite(guess))
        guess = this->Median(); /// this shouldn't be nan or inf
    RealType root = 0;
    math::findMin<RealType>([this](const RealType &x) { return -this->logf(x); }, guess, root);
    return root;
}

template <typename RealType>
long double ContinuousDistribution<RealType>::ExpectedValue(const std::function<double(RealType)> &funPtr,
                                                            RealType minPoint, RealType maxPoint) const
{
    /// attempt to calculate expected value by numerical method
    /// use for distributions w/o explicit formula
    RealType lowerBoundary = minPoint, upperBoundary = maxPoint;
    if (this->isRightBounded())
        lowerBoundary = std::max(minPoint, lowerBoundary);
    if (this->isLeftBounded())
        upperBoundary = std::min(maxPoint, upperBoundary);

    if (lowerBoundary >= upperBoundary)
        return 0.0;

    bool isLeftBoundFinite = std::isfinite(lowerBoundary), isRightBoundFinite = std::isfinite(upperBoundary);

    /// Integrate on finite interval [a, b]
    if (isLeftBoundFinite && isRightBoundFinite) {
        return math::integrate(
            [this, funPtr](double x) {
                double y = funPtr(x);
                return (y == 0.0) ? 0.0 : y * f(x);
            },
            lowerBoundary, upperBoundary);
    }

    /// Integrate on semifinite interval [a, inf)
    if (isLeftBoundFinite) {
        return math::integrate(
            [this, funPtr, lowerBoundary](double x) {
                if (x >= 1.0)
                    return 0.0;
                double denom = 1.0 - x;
                double t = lowerBoundary + x / denom;
                double y = funPtr(t);
                if (y == 0.0)
                    return 0.0;
                y *= this->f(t);
                denom *= denom;
                return y / denom;
            },
            0.0, 1.0);
    }

    /// Integrate on semifinite intervale (-inf, b]
    if (isRightBoundFinite) {
        return math::integrate(
            [this, funPtr, upperBoundary](double x) {
                if (x <= 0.0)
                    return 0.0;
                double t = upperBoundary - (1.0 - x) / x;
                double y = funPtr(t);
                if (y == 0.0)
                    return 0.0;
                y *= this->f(t);
                return y / (x * x);
            },
            0.0, 1.0);
    }

    /// Infinite case
    return math::integrate(
        [this, funPtr](double x) {
            if (std::fabs(x) >= 1.0)
                return 0.0;
            double x2 = x * x;
            double denom = 1.0 - x2;
            double t = x / denom;
            double y = funPtr(t);
            if (y == 0.0)
                return 0.0;
            y *= this->f(t);
            denom *= denom;
            return y * (1.0 + x2) / denom;
        },
        -1.0, 1.0);
}

template <typename RealType>
double ContinuousDistribution<RealType>::Hazard(const RealType &x) const
{
    if (x < this->MinValue())
        return 0.0; /// 0/1
    if (x > this->MaxValue())
        return NAN; /// 0/0
    return this->f(x) / this->S(x);
}

template <typename RealType>
double ContinuousDistribution<RealType>::LikelihoodFunction(std::span<const RealType> sample) const
{
    return std::exp(LogLikelihoodFunction(sample));
}

template <typename RealType>
double ContinuousDistribution<RealType>::LogLikelihoodFunction(std::span<const RealType> sample) const
{
    long double res = 0.0;
    for (const RealType &var : sample)
        res += this->logf(var);
    return res;
}

template <typename RealType>
bool ContinuousDistribution<RealType>::KolmogorovSmirnovTest(std::span<const RealType> orderStatistic,
                                                             double alpha) const
{
    KolmogorovSmirnovRand KSRand;
    double K = KSRand.Quantile1m(alpha);
    size_t n = orderStatistic.size();
    double interval = K / std::sqrt(n);
    double nInv = 1.0 / n;
    double Fn = 0.0;
    for (size_t i = 1; i != n; ++i) {
        RealType x = orderStatistic[i - 1];
        if (x > orderStatistic[i])
            throw std::invalid_argument("Order statistic should be sorted in ascending order");
        double upperBound = Fn + interval;
        Fn = i * nInv;
        double lowerBound = Fn - interval;
        double FReal = this->F(x);
        if (FReal < lowerBound || FReal > upperBound)
            return false;
    }
    double SReal = this->S(orderStatistic[n - 1]);
    return (SReal > interval || SReal < nInv - interval) ? false : true;
}

} // namespace randlib

#include "ContinuousDistribution.h"
#include "KolmogorovSmirnovRand.h"
#include <limits>

template< typename RealType >
void ContinuousDistribution<RealType>::ProbabilityDensityFunction(const std::vector<RealType> &x, std::vector<double> &y) const
{
    for (size_t i = 0; i != x.size(); ++i)
        y[i] = this->f(x[i]);
}

template< typename RealType >
void ContinuousDistribution<RealType>::LogProbabilityDensityFunction(const std::vector<RealType> &x, std::vector<double> &y) const
{
    for (size_t i = 0; i != x.size(); ++i)
        y[i] = this->logf(x[i]);
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl(double p, RealType initValue) const
{
    return searchQuantile(p, false, initValue);
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl(double p) const
{
    return searchQuantile(p, false, std::numeric_limits<RealType>::quiet_NaN());
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl1m(double p, RealType initValue) const
{
    return searchQuantile(p, true, initValue);
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl1m(double p) const
{
    return searchQuantile(p, true, std::numeric_limits<RealType>::quiet_NaN());
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::searchQuantile(double p, bool upperTail, RealType guess) const
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
    auto g = [this, p, upperTail, logScale, logP] (RealType x) {
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
        RealType step = useGuess ? std::max<RealType>(std::fabs(guess) * 1e-3, 1e-3) : 1;
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
        }
        else {
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
    auto split = [anchor, loFinite, hiFinite] (RealType left, RealType right) {
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

    /// Illinois regula falsi, with a bisection step whenever the bracket fails to halve
    int retained = 0;
    bool bisect = false;
    for (int iter = 0; iter < 2000; ++iter) {
        RealType width = b - a, middle;
        if (!bisect && std::isfinite(ga) && std::isfinite(gb) && gb > ga)
            middle = b - gb * (b - a) / (gb - ga);
        else
            middle = split(a, b);
        if (!(middle > a && middle < b))
            middle = split(a, b);
        if (!(middle > a && middle < b))
            break;
        double gm = g(middle);
        if (gm == 0)
            return middle;
        if (gm < 0) {
            a = middle;
            ga = gm;
            if (retained == 1)
                gb *= 0.5;
            retained = 1;
        }
        else {
            b = middle;
            gb = gm;
            if (retained == -1)
                ga *= 0.5;
            retained = -1;
        }
        bisect = (b - a) > 0.5 * width;
    }
    return (std::isfinite(ga) && std::fabs(ga) < std::fabs(gb)) ? a : b;
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::Mode() const
{
    RealType guess = this->Mean(); /// good starting point
    if (!std::isfinite(guess))
        guess = this->Median(); /// this shouldn't be nan or inf
    RealType root = 0;
    RandMath::findMin<RealType>([this] (const RealType &x)
    {
        return -this->logf(x);
    }, guess, root);
    return root;
}

template< typename RealType >
long double ContinuousDistribution<RealType>::ExpectedValue(const std::function<double (RealType)> &funPtr, RealType minPoint, RealType maxPoint) const
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
        return RandMath::integral([this, funPtr] (double x)
        {
            double y = funPtr(x);
            return (y == 0.0) ? 0.0 : y * f(x);
        },
        lowerBoundary, upperBoundary);
    }

    /// Integrate on semifinite interval [a, inf)
    if (isLeftBoundFinite) {
        return RandMath::integral([this, funPtr, lowerBoundary] (double x)
        {
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
        return RandMath::integral([this, funPtr, upperBoundary] (double x)
        {
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
    return RandMath::integral([this, funPtr] (double x)
    {
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
    }, -1.0, 1.0);
}

template< typename RealType >
double ContinuousDistribution<RealType>::Hazard(const RealType &x) const
{
    if (x < this->MinValue())
        return 0.0; /// 0/1
    if (x > this->MaxValue())
        return NAN; /// 0/0
    return this->f(x) / this->S(x);
}

template< typename RealType >
double ContinuousDistribution<RealType>::LikelihoodFunction(const std::vector<RealType> &sample) const
{
    return std::exp(LogLikelihoodFunction(sample));
}

template< typename RealType >
double ContinuousDistribution<RealType>::LogLikelihoodFunction(const std::vector<RealType> &sample) const
{
    long double res = 0.0;
    for (const RealType & var : sample)
        res += this->logf(var);
    return res;
}

template< typename RealType >
bool ContinuousDistribution<RealType>::KolmogorovSmirnovTest(const std::vector<RealType> &orderStatistic, double alpha) const
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

template class ContinuousDistribution<float>;
template class ContinuousDistribution<double>;
template class ContinuousDistribution<long double>;

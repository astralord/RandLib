#include "ContinuousDistribution.h"
#include "KolmogorovSmirnovRand.h"

template< typename RealType >
void ContinuousDistribution<RealType>::ProbabilityDensityFunction(const std::vector<RealType> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = this->f(x[i]);
}

template< typename RealType >
void ContinuousDistribution<RealType>::LogProbabilityDensityFunction(const std::vector<RealType> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = this->logf(x[i]);
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl(double p, RealType initValue) const
{
    static constexpr double SMALL_P = 1e-5;
    if (p < SMALL_P) {
        /// for small p we use logarithmic scale
        double logP = std::log(p);
        if (RandMath::findRoot<RealType>([this, logP] (const RealType & x)
        {
            double logCdf = std::log(this->F(x)), logPdf = this->logf(x);
            double first = logCdf - logP;
            double second = std::exp(logPdf - logCdf);
            return DoublePair(first, second);
        }, initValue))
            throw std::runtime_error("Continuous distribution: failure in numeric procedure");
        return initValue;
    }

    if (this->SupportType() == FINITE_T) {
        if (RandMath::findRoot<RealType>([this, p] (const RealType & x)
        {
            return this->F(x) - p;
        }, this->MinValue(), this->MaxValue(), initValue))
            throw std::runtime_error("Continuous distribution: failure in numeric procedure");
        return initValue;
    }

    if (RandMath::findRoot<RealType>([this, p] (const RealType & x)
    {
        double first = this->F(x) - p;
        double second = this->f(x);
        return DoublePair(first, second);
    }, initValue))
        throw std::runtime_error("Continuous distribution: failure in numeric procedure");
    return initValue;
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl(double p) const
{
    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<RealType> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    if (index == 0.0)
        return this->quantileImpl(p, *std::min_element(sample.begin(), sample.end()));
    std::nth_element(sample.begin(), sample.begin() + index, sample.end());
    return this->quantileImpl(p, sample[index]);
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl1m(double p, RealType initValue) const
{
    static constexpr double SMALL_P = 1e-5;
    if (p < SMALL_P) {
        /// for small p we use logarithmic scale
        double logP = std::log(p);
        if (RandMath::findRoot<RealType>([this, logP] (const RealType & x)
        {
            double logCcdf = std::log(this->S(x)), logPdf = this->logf(x);
            double first = logP - logCcdf;
            double second = std::exp(logPdf - logCcdf);
            return DoublePair(first, second);
        }, initValue))
            throw std::runtime_error("Continuous distribution: failure in numeric procedure");
        return initValue;
    }

    if (this->SupportType() == FINITE_T) {
        if (!RandMath::findRoot<RealType>([this, p] (const RealType & x)
        {
            return this->S(x) - p;
        }, this->MinValue(), this->MaxValue(), initValue))
            throw std::runtime_error("Continuous distribution: failure in numeric procedure");
        return initValue;
    }

    if (!RandMath::findRoot<RealType>([this, p] (const RealType & x)
    {
        double first = p - this->S(x);
        double second = this->f(x);
        return DoublePair(first, second);
    }, initValue))
        throw std::runtime_error("Continuous distribution: failure in numeric procedure");
    return initValue;
}

template< typename RealType >
RealType ContinuousDistribution<RealType>::quantileImpl1m(double p) const
{
    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<RealType> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    if (index == 0.0)
        return this->quantileImpl1m(p, *std::max_element(sample.begin(), sample.end()));
    std::nth_element(sample.begin(), sample.begin() + index, sample.end(), std::greater<>());
    return this->quantileImpl1m(p, sample[index]);
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
    long double res = 1.0;
    for (const RealType & var : sample)
        res *= this->f(var);
    return res;
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
            throw std::invalid_argument("Sample should be sorted in ascending order");
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

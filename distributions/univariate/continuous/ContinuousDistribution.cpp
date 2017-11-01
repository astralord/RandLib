#include "ContinuousDistribution.h"
#include "KolmogorovSmirnovRand.h"

void ContinuousDistribution::ProbabilityDensityFunction(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

void ContinuousDistribution::LogProbabilityDensityFunction(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = x.size();
    if (size > y.size())
        return;
    for (size_t i = 0; i != size; ++i)
        y[i] = logf(x[i]);
}

double ContinuousDistribution::quantileImpl(double p) const
{
    double guess = 0.0;
    SUPPORT_TYPE supp = SupportType();
    if (supp == FINITE_T && p > 1e-5) {
        if (RandMath::findRoot([this, p] (double x)
        {
            return F(x) - p;
        }, MinValue(), MaxValue(), guess))
            return guess;
        return NAN;
    }

    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<double> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    /// if p is too small
    if (index == 0) {
        guess = *std::min_element(sample.begin(), sample.end());
        double logP = std::log(p);
        if (RandMath::findRoot([this, logP] (double x)
        {
            double logCdf = std::log(F(x)), logPdf = logf(x);
            double first = logCdf - logP;
            double second = std::exp(logPdf - logCdf);
            return DoublePair(first, second);
        }, guess))
            return guess;
        return NAN;
    }

    std::nth_element(sample.begin(), sample.begin() + index, sample.end());
    guess = sample[index];

    if (RandMath::findRoot([this, p] (double x)
    {
        double first = F(x) - p;
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double ContinuousDistribution::quantileImpl1m(double p) const
{
    double guess = 0.0;
    SUPPORT_TYPE supp = SupportType();
    /// we use this method only for sufficient large p
    /// in order to avoid underflow
    if (supp == FINITE_T && p > 1e-5) {
        if (RandMath::findRoot([this, p] (double x)
        {
            return S(x) - p;
        }, MinValue(), MaxValue(), guess))
            return guess;
        return NAN;
    }

    /// We use quantile from sample as an initial guess
    static constexpr int SAMPLE_SIZE = 128;
    static std::vector<double> sample(SAMPLE_SIZE);
    this->Sample(sample);
    int index = p * SAMPLE_SIZE;
    /// if p is too small
    if (index == 0) {
        guess = *std::max_element(sample.begin(), sample.end());
        double logP = std::log(p);
        if (RandMath::findRoot([this, logP] (double x)
        {
            double logCcdf = std::log(S(x)), logPdf = logf(x);
            double first = logP - logCcdf;
            double second = std::exp(logPdf - logCcdf);
            return DoublePair(first, second);
        }, guess))
            return guess;
        return NAN;
    }

    /// Partially sort in desceding order
    std::nth_element(sample.begin(), sample.begin() + index, sample.end(), std::greater<>());
    guess = sample[index];

    if (RandMath::findRoot([this, p] (double x)
    {
        double first = p - S(x);
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double ContinuousDistribution::Mode() const
{
    double guess = Mean(); /// good starting point
    if (!std::isfinite(guess))
        guess = Median(); /// this shouldn't be nan or inf
    double root = 0;
    RandMath::findMin([this] (double x)
    {
        return -logf(x);
    }, guess, root);
    return root;
}

double ContinuousDistribution::ExpectedValue(const std::function<double (double)> &funPtr, double minPoint, double maxPoint) const
{
    /// attempt to calculate expected value by numerical method
    /// use for distributions w/o explicit formula
    double lowerBoundary = minPoint, upperBoundary = maxPoint;
    if (isRightBounded())
        lowerBoundary = std::max(minPoint, lowerBoundary);
    if (isLeftBounded())
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
            y *= f(t);
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
            y *= f(t);
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
        y *= f(t);
        denom *= denom;
        return y * (1.0 + x2) / denom;
    }, -1.0, 1.0);
}

double ContinuousDistribution::Hazard(double x) const
{
    if (x < MinValue())
        return 0.0; /// 0/1
    if (x > MaxValue())
        return NAN; /// 0/0
    return f(x) / S(x);
}

double ContinuousDistribution::LikelihoodFunction(const std::vector<double> &sample) const
{
    double res = 1.0;
    for (const double & var : sample)
        res *= f(var);
    return res;
}

double ContinuousDistribution::LogLikelihoodFunction(const std::vector<double> &sample) const
{
    double res = 0.0;
    for (const double & var : sample)
        res += logf(var);
    return res;
}

bool ContinuousDistribution::KolmogorovSmirnovTest(const std::vector<double> &orderStatistic, double alpha) const
{
    KolmogorovSmirnovRand KSRand;
    double K = KSRand.Quantile1m(alpha);
    size_t n = orderStatistic.size();
    double interval = K / std::sqrt(n);
    double nInv = 1.0 / n;
    double Fn = 0.0;
    for (size_t i = 1; i != n; ++i) {
        double x = orderStatistic[i - 1];
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
